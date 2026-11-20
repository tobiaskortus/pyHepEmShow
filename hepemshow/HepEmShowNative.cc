#include "SimulationRun.hh"

#include <dsbo/hdf5.hpp>

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace io = dsbo::io;

namespace {

struct Options {
  std::string input;
  std::string data_file;
  std::string particle = "e-";
  int num_events = 1000;
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
    else if (argument == "--seed") options.seed = std::stod(next_value(index, argc, argv, argument));
    else if (argument == "--particle") options.particle = next_value(index, argc, argv, argument);
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
  return options;
}

struct Design {
  std::vector<double> absorber;
  std::vector<double> gap;
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
  // The simulation consumes doubles; storage reads still preserve file dtypes.
  for (auto* array : {&absorber, &gap})
    if (array->type != io::dtype<double>()) *array = io::convert(array->view(), io::dtype<double>());
  const auto* absorber_data = absorber.data<double>();
  const auto* gap_data = gap.data<double>();
  return {
    std::vector<double>(absorber_data, absorber_data + absorber.size()),
    std::vector<double>(gap_data, gap_data + gap.size()),
  };
}

void require_option(const std::string& value, const std::string& name) {
  if (value.empty()) throw std::runtime_error("missing required option " + name);
}

void validate_settings(const Options& options) {
  require_option(options.data_file, "--data-file");
  if (options.num_events < 0)
    throw std::runtime_error("--num-events must be non-negative");
  if (options.gy <= 0 || options.gz <= 0)
    throw std::runtime_error("--gy and --gz must be positive");
}

int run_forward(const Options& options) {
  validate_settings(options);
  require_option(options.input, "--input");
  require_option(options.edep_output, "--edep-output");
  require_option(options.energy_output, "--energy-output");
  require_option(options.labels_output, "--labels-output");

  const Design design = read_design(options.input);
  const int num_layers = static_cast<int>(design.absorber.size());

  const SimulationRunResult result = simulate(
    options.data_file,
    design.absorber,
    design.gap,
    options.transverse_size,
    options.particle,
    options.energy_min,
    options.energy_max,
    options.num_events,
    options.seed,
    options.verbosity,
    options.gy,
    options.gz,
    options.gaussian_sigma,
    {},
    {},
    {},
    options.absorber_passive);

  io::write_primal(options.energy_output, result.particle_energies, {options.num_events});
  io::write_primal(options.edep_output, result.edep_per_cell,
                  {options.num_events, num_layers, options.gy, options.gz});

  const std::vector<float> labels(static_cast<std::size_t>(options.num_events), 0.0f);
  io::write_primal(options.labels_output, labels);

  if (!options.design_output.empty()) {
    std::vector<double> repeated_design;
    repeated_design.reserve(static_cast<std::size_t>(options.num_events) * 2 * design.absorber.size());
    for (int event = 0; event < options.num_events; ++event) {
      repeated_design.insert(repeated_design.end(), design.absorber.begin(), design.absorber.end());
      repeated_design.insert(repeated_design.end(), design.gap.begin(), design.gap.end());
    }

    io::write_primal(options.design_output, repeated_design, {options.num_events, 2 * num_layers});
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
  const io::Array adjoint = io::load_backward_gradient(options.edep_output, options.edep_adjoint);

  if (adjoint.shape != primal_shape)
    throw std::runtime_error("edep adjoint shape does not match the edep primal shape");

  const std::vector<double> bar_edep(adjoint.data<double>(), adjoint.data<double>() + adjoint.size());

  const SimulationRunResult result = simulate(
    options.data_file,
    design.absorber,
    design.gap,
    options.transverse_size,
    options.particle,
    options.energy_min,
    options.energy_max,
    options.num_events,
    options.seed,
    options.verbosity,
    options.gy,
    options.gz,
    options.gaussian_sigma,
    bar_edep,
    {},
    {},
    options.absorber_passive);

  io::write_gradient(options.input_adjoint, {
    io::field("absorber", result.bar_thickness_absorber, {num_layers}),
    io::field("gap", result.bar_thickness_gap, {num_layers}),
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
    << "          --particle NAME --energy-min VALUE --energy-max VALUE\n"
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
