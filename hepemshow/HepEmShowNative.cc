#include "SimulationRun.hh"

#include <dsbo/hdf5.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace io = dsbo::io;

namespace {

struct Options {
  std::string input;
  std::string data_file;
  std::vector<std::string> particles;
  int num_events = 1000;
  int num_sims_per_particle = 1;
  int num_threads = 0;
  double seed = 1234.0;
  double energy_min = 10000.0;
  double energy_max = 10000.0;
  double transverse_size = 400.0;
  int gy = 1;
  int gz = 1;
  double gaussian_sigma = -1.0;
  int verbosity = 1;
  bool absorber_passive = false;
  std::string edep_output;
  std::string energy_output;
  std::string labels_output;
  std::string design_output;
  std::string edep_adjoint;
  std::string input_adjoint;
};

std::string next_value(int& index, int argc, char** argv, const std::string& option) {
  if (index + 1 >= argc) throw std::runtime_error("missing value for " + option);
  return argv[++index];
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 2; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--input") options.input = next_value(index, argc, argv, argument);
    else if (argument == "--data-file") options.data_file = next_value(index, argc, argv, argument);
    else if (argument == "--num-events") options.num_events = std::stoi(next_value(index, argc, argv, argument));
    else if (argument == "--num-sims-per-particle")
      options.num_sims_per_particle = std::stoi(next_value(index, argc, argv, argument));
    else if (argument == "--num-threads") {
      options.num_threads = std::stoi(next_value(index, argc, argv, argument));
      if (options.num_threads < 1) throw std::runtime_error("--num-threads must be positive");
    }
    else if (argument == "--seed") options.seed = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--particle") options.particles.push_back(next_value(index, argc, argv, argument));
    else if (argument == "--energy-min") options.energy_min = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--energy-max") options.energy_max = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--transverse-size") options.transverse_size = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--gy") options.gy = std::stoi(next_value(index, argc, argv, argument));
    else if (argument == "--gz") options.gz = std::stoi(next_value(index, argc, argv, argument));
    else if (argument == "--gaussian-sigma") options.gaussian_sigma = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--verbosity") options.verbosity = std::stoi(next_value(index, argc, argv, argument));
    else if (argument == "--absorber-passive") options.absorber_passive = true;
    else if (argument == "--edep-output") options.edep_output = next_value(index, argc, argv, argument);
    else if (argument == "--energy-output") options.energy_output = next_value(index, argc, argv, argument);
    else if (argument == "--labels-output") options.labels_output = next_value(index, argc, argv, argument);
    else if (argument == "--design-output") options.design_output = next_value(index, argc, argv, argument);
    else if (argument == "--edep-adjoint") options.edep_adjoint = next_value(index, argc, argv, argument);
    else if (argument == "--input-adjoint") options.input_adjoint = next_value(index, argc, argv, argument);
    else throw std::runtime_error("unknown option: " + argument);
  }
  if (options.particles.empty()) options.particles.push_back("e-");
  return options;
}

struct Design {
  std::vector<double> absorber;
  std::vector<double> gap;
  io::DType output_dtype;
};

Design read_design(const std::string& path) {
  auto values = io::open_dict(path, std::vector<std::string>{"absorber", "gap"});
  auto& absorber = io::at(values, "absorber");
  auto& gap = io::at(values, "gap");
  if (absorber.shape.size() != 1 || gap.shape.size() != 1)
    throw std::runtime_error("input design fields 'absorber' and 'gap' must be one-dimensional");
  if (absorber.shape[0] != gap.shape[0])
    throw std::runtime_error("input design fields 'absorber' and 'gap' must have the same length");
  if (absorber.shape[0] == 0)
    throw std::runtime_error("input design must contain at least one layer");
  // Transport uses doubles, but the Python backend stores results in the absorber dtype.
  const auto output_dtype = absorber.type;
  for (auto* array : {&absorber, &gap})
    if (array->type != io::dtype<double>()) *array = io::convert(array->view(), io::dtype<double>());
  const auto* absorber_data = absorber.data<double>();
  const auto* gap_data = gap.data<double>();
  return {
    std::vector<double>(absorber_data, absorber_data + absorber.size()),
    std::vector<double>(gap_data, gap_data + gap.size()),
    output_dtype,
  };
}

void require_option(const std::string& value, const std::string& name) {
  if (value.empty()) throw std::runtime_error("missing required option " + name);
}

void validate_settings(const Options& options) {
  require_option(options.data_file, "--data-file");
  if (options.num_events < 0)
    throw std::runtime_error("--num-events must be non-negative");
  if (options.num_sims_per_particle < 1)
    throw std::runtime_error("--num-sims-per-particle must be positive");
  if (options.gy <= 0 || options.gz <= 0)
    throw std::runtime_error("--gy and --gz must be positive");
}

struct EventStream {
  std::size_t start;
  int count;
  int seed;
  std::size_t particle_index;
};

std::size_t total_events(const Options& options) {
  return io::elements({options.num_events, options.num_sims_per_particle, options.particles.size()});
}

std::vector<EventStream> make_streams(const Options& options) {
  // URandom accepts a signed 32-bit seed. Preserve valid legacy seeds and wrap larger replica seeds.
  constexpr double seed_period = 4294967296.0;
  double seed = std::fmod(std::trunc(options.seed), seed_period);
  if (seed >= seed_period / 2) seed -= seed_period;
  if (seed < -seed_period / 2) seed += seed_period;
  const int first_seed = static_cast<int>(seed);
  if (options.num_events == 0) return {{0, 0, first_seed, 0}};

  // The stream boundaries and seeds must not depend on scheduling or thread count.
  const int events_per_stream = options.num_threads == 0 ? options.num_events : 16;
  std::mt19937 seeds(static_cast<std::uint32_t>(first_seed));
  std::vector<EventStream> streams;
  std::size_t start = 0;
  for (std::size_t particle = 0; particle < options.particles.size(); ++particle) {
    for (int sim = 0; sim < options.num_sims_per_particle; ++sim) {
      for (int remaining = options.num_events; remaining > 0;) {
        const int count = std::min(events_per_stream, remaining);
        streams.push_back({start, count, start == 0 ? first_seed : static_cast<int>(seeds() & 0x7fffffffU), particle});
        start += count;
        remaining -= count;
      }
    }
  }
  return streams;
}

template <class Run>
void run_streams(const Options& options, const std::vector<EventStream>& streams, Run run) {
  // Physics tables are passive, read-only inputs. Each simulate() call owns its mutable state.
  const auto state = load_simulation_state(options.data_file);
  const auto workers = std::min(static_cast<std::size_t>(std::max(1, options.num_threads)), streams.size());

  // With one worker this runs sequentially in stream order; with several it is a work queue.
  std::atomic<std::size_t> next{0};
  auto worker = [&] {
    for (auto index = next.fetch_add(1, std::memory_order_relaxed); index < streams.size();
         index = next.fetch_add(1, std::memory_order_relaxed))
      run(index, streams[index], state.get());
  };
  std::vector<std::future<void>> threads;
  threads.reserve(workers - 1);
  for (std::size_t i = 1; i < workers; ++i)
    threads.push_back(std::async(std::launch::async, worker));
  worker();
  for (auto& thread : threads) thread.get();
}

SimulationRunResult simulate_stream(
    const Options& options, const Design& design, const EventStream& stream, const G4HepEmState* state,
    const std::vector<double>& adjoint = {}) {
  return simulate(
    options.data_file, design.absorber, design.gap, options.transverse_size, options.particles[stream.particle_index],
    options.energy_min, options.energy_max, stream.count, stream.seed,
    options.num_threads == 0 ? options.verbosity : 0, options.gy, options.gz, options.gaussian_sigma,
    adjoint, {}, {}, options.absorber_passive, state);
}

void write_output(const std::string& path, io::ArrayView values, io::DType dtype) {
  if (values.type == dtype) io::HDF5Storage::write_primal(path, values);
  else io::HDF5Storage::write_primal(path, io::convert(values, dtype).view());
}

int run_forward(const Options& options) {
  validate_settings(options);
  require_option(options.input, "--input");
  require_option(options.edep_output, "--edep-output");
  require_option(options.energy_output, "--energy-output");
  require_option(options.labels_output, "--labels-output");

  const Design design = read_design(options.input);
  const auto num_layers = design.absorber.size();
  const auto events = total_events(options);

  const std::size_t cells = design.absorber.size() * options.gy * options.gz;
  const io::Shape edep_shape{events, num_layers, options.gy, options.gz};
  std::vector<double> energies(events), edep(io::elements(edep_shape));
  std::vector<float> labels(events);
  const auto streams = make_streams(options);
  run_streams(options, streams, [&](std::size_t, const EventStream& stream, const G4HepEmState* state) {
    const auto result = simulate_stream(options, design, stream, state);
    std::copy(result.particle_energies.begin(), result.particle_energies.end(), energies.begin() + stream.start);
    std::copy(result.edep_per_cell.begin(), result.edep_per_cell.end(), edep.begin() + stream.start * cells);
    std::fill_n(labels.begin() + stream.start, stream.count, static_cast<float>(stream.particle_index));
  });

  write_output(options.energy_output, io::view_of(energies), design.output_dtype);
  write_output(options.edep_output, io::view_of(edep, edep_shape), design.output_dtype);
  io::write_primal(options.labels_output, labels);

  if (!options.design_output.empty()) {
    std::vector<double> repeated_design;
    const io::Shape design_shape{events, 2 * num_layers};
    repeated_design.reserve(io::elements(design_shape));
    for (std::size_t event = 0; event < events; ++event) {
      repeated_design.insert(repeated_design.end(), design.absorber.begin(), design.absorber.end());
      repeated_design.insert(repeated_design.end(), design.gap.begin(), design.gap.end());
    }

    io::write_primal(options.design_output, repeated_design, design_shape);
  }

  return 0;
}

#ifdef CODI_REVERSE

int run_backward(const Options& options) {
  validate_settings(options);
  require_option(options.input, "--input");
  require_option(options.edep_output, "--edep-output");
  require_option(options.edep_adjoint, "--edep-adjoint");
  require_option(options.input_adjoint, "--input-adjoint");

  const Design design = read_design(options.input);
  const int num_layers = static_cast<int>(design.absorber.size());

  const io::Shape primal_shape = io::shape_from_path(options.edep_output);
  io::Array adjoint = io::load_backward_gradient(options.edep_output, options.edep_adjoint);

  if (adjoint.shape != primal_shape)
    throw std::runtime_error("edep adjoint shape does not match the edep primal shape");

  if (primal_shape != io::Shape{total_events(options), num_layers, options.gy, options.gz})
    throw std::runtime_error("edep shape does not match the event count and geometry");
  if (adjoint.type != io::dtype<double>()) adjoint = io::convert(adjoint.view(), io::dtype<double>());

  const std::size_t cells = design.absorber.size() * options.gy * options.gz;
  const auto streams = make_streams(options);
  std::vector<std::vector<double>> absorber_parts(streams.size()), gap_parts(streams.size());
  run_streams(options, streams, [&](std::size_t index, const EventStream& stream, const G4HepEmState* state) {
    std::vector<double> weights;
    if (stream.count > 0) {
      const auto* begin = adjoint.data<double>() + stream.start * cells;
      weights.assign(begin, begin + stream.count * cells);
    }
    auto result = simulate_stream(options, design, stream, state, weights);
    absorber_parts[index] = std::move(result.bar_thickness_absorber);
    gap_parts[index] = std::move(result.bar_thickness_gap);
  });

  std::vector<double> absorber(num_layers), gap(num_layers);
  // Sum in stream order so scheduling does not change floating-point accumulation.
  for (std::size_t index = 0; index < streams.size(); ++index)
    for (int layer = 0; layer < num_layers; ++layer) {
      absorber[layer] += absorber_parts[index][layer];
      gap[layer] += gap_parts[index][layer];
    }

  io::write_gradient(options.input_adjoint, {
    io::field("absorber", absorber, {num_layers}),
    io::field("gap", gap, {num_layers}),
  });

  return 0;
}

#else

int run_backward(const Options&) {
  std::cerr << "hepemshow-native: backward requires a reverse-mode build (CODI_REVERSE)." << std::endl;
  return 1;
}

#endif

void print_usage() {
  std::cerr
    << "usage: hepemshow-native <forward|backward> [options]\n"
    << "  common: --input PATH --data-file PATH --num-events N --seed VALUE\n"
    << "          [--num-threads N]  use deterministic event streams with up to N threads\n"
    << "          [--particle NAME]... [--num-sims-per-particle N]\n"
    << "          --energy-min VALUE --energy-max VALUE\n"
    << "          --transverse-size VALUE --gy N --gz N --gaussian-sigma VALUE\n"
    << "          --absorber-passive --verbosity N\n"
    << "  forward: --edep-output PATH --energy-output PATH --labels-output PATH\n"
    << "           [--design-output PATH]\n"
    << "  backward: --edep-output PATH --edep-adjoint PATH --input-adjoint PATH\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    io::check_environment();
    if (argc < 2) {
      print_usage();
      return 2;
    }

    const std::string mode = argv[1];
    const Options options = parse_options(argc, argv);

    if (mode == "forward") return run_forward(options);
    if (mode == "backward") return run_backward(options);

    std::cerr << "hepemshow-native: unknown mode '" << mode << "'" << std::endl;
    print_usage();
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "hepemshow-native: " << error.what() << std::endl;
    return 1;
  }
}
