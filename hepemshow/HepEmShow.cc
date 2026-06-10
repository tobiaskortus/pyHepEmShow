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
    const std::vector<double>& bar_edep               = {},
    // CODI_FORWARD tangent seeds (per-layer or empty for zero; ignored in non-AD builds)
    const std::vector<double>& absorber_thickness_dot = {},
    const std::vector<double>& gap_thickness_dot      = {})
{
  // Derive number of layers from the thickness vectors
  if (absorber_thickness.size() != gap_thickness.size())
    throw std::runtime_error(
      "run_simulation: absorber_thickness and gap_thickness must have the same length "
      "(got " + std::to_string(absorber_thickness.size()) + " vs "
              + std::to_string(gap_thickness.size()) + ")");
              
  const int num_layers = (int)absorber_thickness.size();

  if (!absorber_thickness_dot.empty() && (int)absorber_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: absorber_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(absorber_thickness_dot.size()));

  if (!gap_thickness_dot.empty() && (int)gap_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: gap_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(gap_thickness_dot.size()));

#ifdef CODI_REVERSE
  const std::size_t num_layers_sz = static_cast<std::size_t>(num_layers);
  const std::size_t num_events_sz = static_cast<std::size_t>(num_events);
  const std::size_t event_layer_size = num_events_sz * num_layers_sz;
  if (!bar_edep.empty() && bar_edep.size() != num_layers_sz && bar_edep.size() != event_layer_size) {
    throw std::runtime_error(
      "run_simulation: bar_edep must be empty, per-layer (" + std::to_string(num_layers_sz) +
      "), or per-event-per-layer (" + std::to_string(event_layer_size) + "), got " +
      std::to_string(bar_edep.size()));
  }
#endif
  std::vector<double> res_edep_per_layer;
  std::vector<double> res_particle_energy;
#if CODI_FORWARD
  std::vector<double> res_edep_per_layer_dot;
#endif
#ifdef CODI_REVERSE
  const std::vector<double> bar_edep_local = bar_edep;
  std::vector<double> res_bar_thickness_absorber;
  std::vector<double> res_bar_thickness_gap;
  double              res_bar_particle_energy_mean = 0.0;
  double              res_bar_particle_energy_var  = 0.0;
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
    thePrimaryGenerator.SetKinEnergy(particle_energy_min);

    Results theResult;
    theResult.fEdepPerLayer_CurrentEvent.ReSet("hist_Edep_PerLayer_CurrentEvent", 0, num_layers, num_layers);
    theResult.fEdepPerLayer.reserve(num_events * num_layers);
    theResult.fParticleEnergy.reserve(num_events);
#if CODI_FORWARD
    theResult.fEdepPerLayerD.reserve(num_events * num_layers);
#endif
#ifdef CODI_REVERSE
    theResult.barEdep = bar_edep_local.empty()
      ? std::vector<double>(static_cast<std::size_t>(num_events) * static_cast<std::size_t>(num_layers), 0.0)
      : bar_edep_local;
    theResult.barThicknessAbsorber.resize(num_layers);
    theResult.barThicknessGap.resize(num_layers);
#endif
    EventLoop::ProcessEvents(*theTLData, *theState, thePrimaryGenerator, theGeometry,
                             theResult, num_events, verbosity,
                             particle_energy_min, particle_energy_max, theURnd);

    delete theRandomEngine;
    delete theURnd;
    delete theTLData;
    FreeG4HepEmData(theState->fData);
    delete theState;

    // Move data before GIL is re-acquired
    res_edep_per_layer   = std::move(theResult.fEdepPerLayer);
    res_particle_energy  = std::move(theResult.fParticleEnergy);
#if CODI_FORWARD
    res_edep_per_layer_dot = std::move(theResult.fEdepPerLayerD);
#endif
#ifdef CODI_REVERSE
    res_bar_thickness_absorber.resize(num_layers);
    res_bar_thickness_gap.resize(num_layers);
    for (int i = 0; i < num_layers; ++i) {
      res_bar_thickness_absorber[i] = theResult.barThicknessAbsorber[i].getMean();
      res_bar_thickness_gap[i]      = theResult.barThicknessGap[i].getMean();
    }
    res_bar_particle_energy_mean = theResult.barParticleEnergy.getMean();
    res_bar_particle_energy_var  = theResult.barParticleEnergy.getVar();
#endif
  } // GIL re-acquired here

  auto to_numpy2d = [num_events, num_layers](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[2] = {(size_t)num_events, (size_t)num_layers};
    return nb::ndarray<nb::numpy, double, nb::ndim<2>>(heap->data(), 2, shape, owner);
  };

  auto to_numpy1d = [num_events](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[1] = {(size_t)num_events};
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(heap->data(), 1, shape, owner);
  };

  nb::dict result;
  result["edep_per_layer"]    = to_numpy2d(std::move(res_edep_per_layer));
  result["particle_energies"] = to_numpy1d(std::move(res_particle_energy));
#if CODI_FORWARD
  result["edep_per_layer_dot"] = to_numpy2d(std::move(res_edep_per_layer_dot));
#endif
#ifdef CODI_REVERSE
  auto to_numpy1d_vec = [num_layers](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
    size_t shape[1] = {(size_t)num_layers};
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(heap->data(), 1, shape, owner);
  };
  result["bar_thickness_absorber"] = to_numpy1d_vec(std::move(res_bar_thickness_absorber));
  result["bar_thickness_gap"]      = to_numpy1d_vec(std::move(res_bar_thickness_gap));
  nb::dict bar_pe;
  bar_pe["mean"] = res_bar_particle_energy_mean;
  bar_pe["var"]  = res_bar_particle_energy_var;
  result["bar_particle_energy"] = bar_pe;
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
      nb::arg("particle_energy_min") = 10000.0,
      nb::arg("particle_energy_max") = 10000.0,
      nb::arg("num_events")         = 1000,
      nb::arg("random_seed")        = 1234.0,
      nb::arg("verbosity")          = 1,
      nb::arg("bar_edep")           = std::vector<double>{},
      nb::arg("absorber_thickness_dot") = std::vector<double>{},
      nb::arg("gap_thickness_dot")      = std::vector<double>{});
  };

  m.attr("compiled_mode") = PYHEPEMSHOW_MODE_NAME;
  bind_run(m);

  nb::module_ mode_ns = m.def_submodule(PYHEPEMSHOW_MODE_NAME);
  bind_run(mode_ns);
}
