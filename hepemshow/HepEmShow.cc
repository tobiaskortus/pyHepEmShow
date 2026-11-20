#include "SimulationRun.hh"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <string>
#include <vector>

namespace nb = nanobind;

#ifndef PYHEPEMSHOW_MODULE_NAME
#define PYHEPEMSHOW_MODULE_NAME _core
#endif


nb::dict run_simulation(
    const std::string&         data_file,
    const std::vector<double>& absorber_thickness     = {},
    const std::vector<double>& gap_thickness          = {},
    double                     transverse_size        = 400.0,
    const std::string&         particle               = "e-",
    double                     particle_energy_min    = 10000.0,
    double                     particle_energy_max    = 10000.0,
    int                        num_events             = 1000,
    double                     random_seed            = 1234.0,
    int                        verbosity              = 1,
    int                        gy          = 1,
    int                        gz          = 1,
    double                     gaussian_sigma         = -1.0,
    const std::vector<double>& bar_edep               = {},
    // CODI_FORWARD tangent seeds (per-layer or empty for zero; ignored in non-AD builds)
    const std::vector<double>& absorber_thickness_dot = {},
    const std::vector<double>& gap_thickness_dot      = {},
    bool                       absorber_passive       = false)
{
  const int num_layers = (int)absorber_thickness.size();
  SimulationRunResult result;
  {
    nb::gil_scoped_release release;
    result = simulate(
      data_file,
      absorber_thickness,
      gap_thickness,
      transverse_size,
      particle,
      particle_energy_min,
      particle_energy_max,
      num_events,
      random_seed,
      verbosity,
      gy,
      gz,
      gaussian_sigma,
      bar_edep,
      absorber_thickness_dot,
      gap_thickness_dot,
      absorber_passive);
  }

  auto to_numpy4d = [num_events, num_layers, gy, gz](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[4] = {(size_t)num_events, (size_t)num_layers, (size_t)gy, (size_t)gz};
    return nb::ndarray<nb::numpy, double, nb::ndim<4>>(heap->data(), 4, shape, owner);
  };

  auto to_numpy1d = [num_events](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[1] = {(size_t)num_events};
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(heap->data(), 1, shape, owner);
  };

  nb::dict out;
  out["edep_per_cell"]    = to_numpy4d(std::move(result.edep_per_cell));
  out["particle_energies"] = to_numpy1d(std::move(result.particle_energies));
#if CODI_FORWARD
  out["edep_per_cell_dot"] = to_numpy4d(std::move(result.edep_per_cell_dot));
#endif
#ifdef CODI_REVERSE
  auto to_numpy1d_vec = [](std::vector<double> v) {
    const size_t num_layers = v.size();
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
    size_t shape[1] = {num_layers};
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(heap->data(), 1, shape, owner);
  };
  out["bar_thickness_absorber"] = to_numpy1d_vec(std::move(result.bar_thickness_absorber));
  out["bar_thickness_gap"]      = to_numpy1d_vec(std::move(result.bar_thickness_gap));
  nb::dict bar_pe;
  bar_pe["mean"] = result.bar_particle_energy_mean;
  bar_pe["var"]  = result.bar_particle_energy_var;
  out["bar_particle_energy"] = bar_pe;
#endif
  return out;
}


NB_MODULE(PYHEPEMSHOW_MODULE_NAME, m) {
  auto bind_run = [](nb::module_& target) {
    target.def("run_simulation", &run_simulation,
      nb::arg("data_file"),
      nb::arg("absorber_thickness") = std::vector<double>(50, 2.3),
      nb::arg("gap_thickness")      = std::vector<double>(50, 5.7),
      nb::arg("transverse_size")    = 400.0,
      nb::arg("particle")           = "e-",
      nb::arg("particle_energy_min") = 10000.0,
      nb::arg("particle_energy_max") = 10000.0,
      nb::arg("num_events")         = 1000,
      nb::arg("random_seed")        = 1234.0,
      nb::arg("verbosity")          = 1,
      nb::arg("gy")      = 1,
      nb::arg("gz")      = 1,
      nb::arg("gaussian_sigma")     = -1.0,
      nb::arg("bar_edep")           = std::vector<double>{},
      nb::arg("absorber_thickness_dot") = std::vector<double>{},
      nb::arg("gap_thickness_dot")      = std::vector<double>{},
      nb::arg("absorber_passive")       = false);
  };

  m.attr("compiled_mode") = PYHEPEMSHOW_MODE_NAME;
  bind_run(m);

  nb::module_ mode_ns = m.def_submodule(PYHEPEMSHOW_MODE_NAME);
  bind_run(mode_ns);
}
