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

namespace nb = nanobind;


nb::dict run_simulation(
    const std::string&      data_file,
    int                     num_layers         = 50,
    double                  absorber_thickness = 2.3,
    double                  gap_thickness      = 5.7,
    double                  transverse_size    = 400.0,
    const std::string&      particle           = "e-",
    double                  particle_energy    = 10000.0,
    int                     num_events         = 1000,
    double                  random_seed        = 1234.0,
    int                     verbosity          = 1,
    const std::vector<double>& bar_edep         = {},
    // CODI_FORWARD tangent seeds (ignored in non-AD builds)
    double                  absorber_thickness_dot = 0.0,
    double                  gap_thickness_dot      = 0.0,
    double                  particle_energy_dot    = 0.0)
{
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
    G4double absThick = absorber_thickness, gapThick = gap_thickness;
#if CODI_FORWARD
    SET_DOTVALUE(absThick, absorber_thickness_dot);
    SET_DOTVALUE(gapThick, gap_thickness_dot);
#endif
    theGeometry.SetAbsThick(absThick);
    theGeometry.SetGapThick(gapThick);
  }

  PrimaryGenerator thePrimaryGenerator;
  thePrimaryGenerator.SetCharge((particle == "e-") ? -1.0 : (particle == "gamma" ? 0.0 : +1.0));
  {
    G4double pEnergy = particle_energy;
#if CODI_FORWARD
    SET_DOTVALUE(pEnergy, particle_energy_dot);
#endif
    thePrimaryGenerator.SetKinEnergy(pEnergy);
  }
  thePrimaryGenerator.SetPosition(theGeometry.GetPrimaryXposition(), 0.0, 0.0);
  thePrimaryGenerator.SetDirection(1.0, 0.0, 0.0);

  Results theResult;
  theResult.fEdepPerLayer.ReSet("hist_Edep_PerLayer", 0, num_layers, num_layers);
  theResult.fEdepPerLayer_CurrentEvent.ReSet("hist_Edep_PerLayer_CurrentEvent", 0, num_layers, num_layers);
  theResult.fEdepPerLayer_Acc.resize(50);
  theResult.fEdepPerLayer_Raw.reserve(num_events * 50);
#if CODI_FORWARD
  theResult.fEdepPerLayer_AccD.resize(50);
  theResult.fEdepPerLayer_RawD.reserve(num_events * 50);
#endif
#ifdef CODI_REVERSE
  theResult.barEdep.resize(50, 0.0);
  for (std::size_t i = 0; i < bar_edep.size() && i < 50; ++i)
    theResult.barEdep[i] = bar_edep[i];
#endif
  theResult.fGammaTrackLenghtPerLayer.ReSet("hist_GamTrackL_PerLayer", 0, num_layers, num_layers);
  theResult.fElPosTrackLenghtPerLayer.ReSet("hist_ElPosTrackL_PerLayer", 0, num_layers, num_layers);

  EventLoop::ProcessEvents(*theTLData, *theState, thePrimaryGenerator, theGeometry,
                           theResult, num_events, verbosity);

  delete theRandomEngine;
  delete theURnd;
  delete theTLData;
  FreeG4HepEmData(theState->fData);
  delete theState;

  auto to_numpy2d = [num_events](std::vector<double> v) {
    auto* heap = new std::vector<double>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept {
      delete static_cast<std::vector<double>*>(p);
    });
    size_t shape[2] = {(size_t)num_events, 50u};
    return nb::ndarray<nb::numpy, double, nb::ndim<2>>(heap->data(), 2, shape, owner);
  };

  nb::dict result;
  result["edep_per_layer"] = to_numpy2d(std::move(theResult.fEdepPerLayer_Raw));
#if CODI_FORWARD
  result["edep_per_layer_dot"] = to_numpy2d(std::move(theResult.fEdepPerLayer_RawD));
#endif
#ifdef CODI_REVERSE
  auto bar = [](const Accumulator<double>& a) -> nb::dict {
    nb::dict d; d["mean"] = a.getMean(); d["var"] = a.getVar(); return d;
  };
  result["bar_thickness_absorber"] = bar(theResult.barThicknessAbsorber);
  result["bar_thickness_gap"]      = bar(theResult.barThicknessGap);
  result["bar_particle_energy"]    = bar(theResult.barParticleEnergy);
#endif
  return result;
}


NB_MODULE(pyhepemshow, m) {
  m.def("run_simulation", &run_simulation,
    nb::arg("data_file"),
    nb::arg("num_layers")         = 50,
    nb::arg("absorber_thickness") = 2.3,
    nb::arg("gap_thickness")      = 5.7,
    nb::arg("transverse_size")    = 400.0,
    nb::arg("particle")           = "e-",
    nb::arg("particle_energy")    = 10000.0,
    nb::arg("num_events")         = 1000,
    nb::arg("random_seed")        = 1234.0,
    nb::arg("verbosity")          = 1,
    nb::arg("bar_edep")           = std::vector<double>{},
    nb::arg("absorber_thickness_dot") = 0.0,
    nb::arg("gap_thickness_dot")      = 0.0,
    nb::arg("particle_energy_dot")    = 0.0);
}
