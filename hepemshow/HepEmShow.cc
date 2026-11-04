#include "ad_type.h"

// G4HepEm includes
#include "G4HepEmState.hh"
#include "G4HepEmData.hh"
#include "G4HepEmDataJsonIO.hh"
#include "G4HepEmTLData.hh"
#include "G4HepEmRandomEngine.hh"

// Local simulation includes
#include "URandom.hh"
#include "Geometry.hh"
#include "PrimaryGenerator.hh"
#include "Results.hh"
#include "EventLoop.hh"

// nanobind
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <fstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <mutex>

namespace nb = nanobind;

using ConstNdArray4D = nb::ndarray<nb::numpy, const double, nb::ndim<4>, nb::c_contig>;

#ifndef PYHEPEMSHOW_MODULE_NAME
#define PYHEPEMSHOW_MODULE_NAME _core
#endif


nb::dict run_simulation(
    const std::string&         data_file,
    const std::vector<double>& absorber_thickness     = {},
    const std::vector<double>& gap_thickness          = {},
    double                     transverse_size        = 400.0,
    const std::string&         particle               = "e-",
    const std::vector<double>& particle_energies      = {},
    double                     random_seed            = 1234.0,
    int                        verbosity              = 1,
    int                        gy          = 1,
    int                        gz          = 1,
    double                     gaussian_sigma         = -1.0,
    const std::string&         kernel      = "trilinear",
    ConstNdArray4D             bar_edep               = {},
    // CODI_FORWARD tangent seeds (per-layer or empty for zero; ignored in non-AD builds)
    const std::vector<double>& absorber_thickness_dot = {},
    const std::vector<double>& gap_thickness_dot      = {},
    const std::vector<double>& particle_energies_dot  = {},
    const std::vector<std::vector<double>>& absorber_thickness_dots = {},
    const std::vector<std::vector<double>>& gap_thickness_dots      = {},
    bool                       absorber_passive       = false)
{
  // Derive number of layers from the thickness vectors
  if (absorber_thickness.size() != gap_thickness.size())
    throw std::runtime_error(
      "run_simulation: absorber_thickness and gap_thickness must have the same length "
      "(got " + std::to_string(absorber_thickness.size()) + " vs "
              + std::to_string(gap_thickness.size()) + ")");
              
  const int num_layers = (int)absorber_thickness.size();
  const int num_events = (int)particle_energies.size();

  if (!absorber_thickness_dot.empty() && (int)absorber_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: absorber_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(absorber_thickness_dot.size()));

  if (!gap_thickness_dot.empty() && (int)gap_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: gap_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(gap_thickness_dot.size()));

  if (!particle_energies_dot.empty() && (int)particle_energies_dot.size() != num_events)
    throw std::runtime_error(
      "run_simulation: particle_energies_dot must be empty or match num_events (" +
      std::to_string(num_events) + "), got " +
      std::to_string(particle_energies_dot.size()));

  auto validate_per_event_layer_dots = [&](const std::vector<std::vector<double>>& dots, const char* name) {
    if (dots.empty())
      return;
    if ((int)dots.size() != num_events)
      throw std::runtime_error(
        std::string("run_simulation: ") + name + " must be empty or have num_events (" +
        std::to_string(num_events) + ") rows, got " + std::to_string(dots.size()));
    for (std::size_t i = 0; i < dots.size(); ++i) {
      if ((int)dots[i].size() != num_layers)
        throw std::runtime_error(
          std::string("run_simulation: ") + name + " row " + std::to_string(i) +
          " must have num_layers (" + std::to_string(num_layers) + ") columns, got " +
          std::to_string(dots[i].size()));
    }
  };

  validate_per_event_layer_dots(absorber_thickness_dots, "absorber_thickness_dots");
  validate_per_event_layer_dots(gap_thickness_dots, "gap_thickness_dots");

  const std::size_t total_cells_per_layer = static_cast<std::size_t>(gy) * static_cast<std::size_t>(gz);
  const std::size_t total_cells_per_event = static_cast<std::size_t>(num_layers) * total_cells_per_layer;

#ifdef CODI_REVERSE
  const std::size_t num_events_sz = static_cast<std::size_t>(num_events);
  const std::size_t event_cell_size = num_events_sz * total_cells_per_event;
  std::vector<double> bar_edep_local;
  
  if (bar_edep.is_valid()) {
    if (bar_edep.shape(0) != static_cast<size_t>(num_events) ||
        bar_edep.shape(1) != static_cast<size_t>(num_layers) ||
        bar_edep.shape(2) != static_cast<size_t>(gy) ||
        bar_edep.shape(3) != static_cast<size_t>(gz)) {
      throw std::runtime_error(
        "run_simulation: bar_edep must have shape (" +
        std::to_string(num_events) + ", " +
        std::to_string(num_layers) + ", " +
        std::to_string(gy) + ", " +
        std::to_string(gz) + ")");
    }
    bar_edep_local.assign(bar_edep.data(), bar_edep.data() + bar_edep.size());
  }
#endif
  std::vector<double> res_edep_per_cell;
  std::vector<double> res_particle_energy;
#if CODI_FORWARD
  std::vector<double> res_edep_per_cell_dot;
#endif
#ifdef CODI_REVERSE
  std::vector<double> res_bar_thickness_absorber;
  std::vector<double> res_bar_thickness_gap;
  std::vector<double> res_bar_particle_energies;
#endif

  {
    nb::gil_scoped_release release;

    std::ifstream jsonIS(data_file.c_str());
    if (!jsonIS.is_open())
      throw std::runtime_error("run_simulation: cannot open data file: " + data_file);
    G4HepEmState* theState = G4HepEmStateFromJson(jsonIS);

    URandom*             theURnd         = new URandom(random_seed);
    G4HepEmRandomEngine* theRandomEngine = new G4HepEmRandomEngine(theURnd);
    G4HepEmTLData*       theTLData       = new G4HepEmTLData();
    theTLData->SetRandomEngine(theRandomEngine);

    Geometry theGeometry;
    theGeometry.SetNumLayers(num_layers);
    theGeometry.SetCaloSizeYZ(transverse_size);
    theGeometry.SetGranularity(gy, gz);
    theGeometry.SetGaussianSigma(gaussian_sigma);
    
    if (kernel == "hard") {
      theGeometry.SetEdepKernel(EdepKernel::Hard);
    } else if (kernel == "gaussian") {
      theGeometry.SetEdepKernel(EdepKernel::Gaussian);
    } else if (kernel == "trilinear") {
      theGeometry.SetEdepKernel(EdepKernel::Trilinear);
    } else {
      throw std::runtime_error(
        "run_simulation: unknown kernel '" + kernel + "' (expected 'hard', 'gaussian' or 'trilinear')");
    }
    
    {
      std::vector<G4double> absThickVec(absorber_thickness.begin(), absorber_thickness.end());
      std::vector<G4double> gapThickVec(gap_thickness.begin(), gap_thickness.end());
#if CODI_FORWARD
      if (!absorber_thickness_dot.empty())
        for (int i = 0; i < num_layers; ++i) SET_DOTVALUE(absThickVec[i], absorber_thickness_dot[i]);
      if (!gap_thickness_dot.empty())
        for (int i = 0; i < num_layers; ++i) SET_DOTVALUE(gapThickVec[i], gap_thickness_dot[i]);
#endif
      theGeometry.SetAbsThickPerLayer(absThickVec);
      theGeometry.SetGapThickPerLayer(gapThickVec);
    }

    PrimaryGenerator thePrimaryGenerator;
    thePrimaryGenerator.SetCharge((particle == "e-") ? -1.0 : (particle == "gamma" ? 0.0 : +1.0));
    thePrimaryGenerator.SetPosition(theGeometry.GetPrimaryXposition(), 0.0, 0.0);
    thePrimaryGenerator.SetDirection(1.0, 0.0, 0.0);

    Results theResult;
    theResult.fAbsorberPassive = absorber_passive;
    theResult.fEdepPerCell_CurrentEvent.assign(num_layers,
      std::vector<std::vector<G4double>>(gy, std::vector<G4double>(gz, 0.0)));
    const std::size_t cellCount = static_cast<std::size_t>(num_layers) * static_cast<std::size_t>(gy) * static_cast<std::size_t>(gz);
    theResult.fEdepPerCell.reserve(num_events * cellCount);
    theResult.fParticleEnergy.reserve(num_events);
#if CODI_FORWARD
    theResult.fEdepPerCellD.reserve(num_events * cellCount);
#endif
#ifdef CODI_REVERSE
    theResult.barEdep = bar_edep_local.empty()
      ? std::vector<double>(event_cell_size, 0.0)
      : bar_edep_local;
    theResult.barThicknessAbsorber.resize(num_layers);
    theResult.barThicknessGap.resize(num_layers);
#endif
    EventLoop::ProcessEvents(*theTLData, *theState, thePrimaryGenerator, theGeometry,
                             theResult, num_events, verbosity,
                             particle_energies, particle_energies_dot,
                             absorber_thickness_dots, gap_thickness_dots,
                             theURnd, static_cast<long>(random_seed));

    delete theRandomEngine;
    delete theURnd;
    delete theTLData;
    FreeG4HepEmData(theState->fData);
    delete theState;

    // Move data before GIL is re-acquired
    res_edep_per_cell   = std::move(theResult.fEdepPerCell);
    res_particle_energy  = std::move(theResult.fParticleEnergy);
#if CODI_FORWARD
    res_edep_per_cell_dot = std::move(theResult.fEdepPerCellD);
#endif
#ifdef CODI_REVERSE
    res_bar_thickness_absorber = std::move(theResult.barThicknessAbsorber);
    res_bar_thickness_gap      = std::move(theResult.barThicknessGap);
    res_bar_particle_energies  = std::move(theResult.barParticleEnergies);
#endif
  } // GIL re-acquired here

  auto to_numpy4d = [num_events, num_layers, gy, gz](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[4] = {(size_t)num_events, (size_t)num_layers, (size_t)gy, (size_t)gz};
    return nb::ndarray<nb::numpy, double, nb::ndim<4>>(heap->data(), 4, shape, owner);
  };

  auto to_numpy1d = [](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[1] = {heap->size()};
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(heap->data(), 1, shape, owner);
  };

  nb::dict result;
  result["edep_per_cell"]    = to_numpy4d(std::move(res_edep_per_cell));
  result["particle_energies"] = to_numpy1d(std::move(res_particle_energy));
#if CODI_FORWARD
  result["edep_per_cell_dot"] = to_numpy4d(std::move(res_edep_per_cell_dot));
#endif
#ifdef CODI_REVERSE
  result["bar_thickness_absorber"] = to_numpy1d(std::move(res_bar_thickness_absorber));
  result["bar_thickness_gap"]      = to_numpy1d(std::move(res_bar_thickness_gap));
  result["bar_particle_energies"]  = to_numpy1d(std::move(res_bar_particle_energies));
#endif
  return result;
}


NB_MODULE(PYHEPEMSHOW_MODULE_NAME, m) {
  auto bind_run = [](nb::module_& target) {
    target.def("run_simulation", &run_simulation,
      nb::arg("data_file"),
      nb::arg("absorber_thickness") = std::vector<double>(50, 2.3),
      nb::arg("gap_thickness")      = std::vector<double>(50, 5.7),
      nb::arg("transverse_size")    = 400.0,
      nb::arg("particle")           = "e-",
      nb::arg("particle_energies"),
      nb::arg("random_seed")        = 1234.0,
      nb::arg("verbosity")          = 1,
      nb::arg("gy")      = 1,
      nb::arg("gz")      = 1,
      nb::arg("gaussian_sigma")     = -1.0,
      nb::arg("kernel")             = "trilinear",
      nb::arg("bar_edep")           = nb::none(),
      nb::arg("absorber_thickness_dot") = std::vector<double>{},
      nb::arg("gap_thickness_dot")      = std::vector<double>{},
      nb::arg("particle_energies_dot")  = std::vector<double>{},
      nb::arg("absorber_thickness_dots") = std::vector<std::vector<double>>{},
      nb::arg("gap_thickness_dots")      = std::vector<std::vector<double>>{},
      nb::arg("absorber_passive")       = false);
  };

  m.attr("compiled_mode") = PYHEPEMSHOW_MODE_NAME;
  bind_run(m);

  nb::module_ mode_ns = m.def_submodule(PYHEPEMSHOW_MODE_NAME);
  bind_run(mode_ns);
}
