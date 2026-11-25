#include "SimulationRun.hh"

#include "G4HepEmState.hh"
#include "G4HepEmData.hh"
#include "G4HepEmDataJsonIO.hh"
#include "G4HepEmParameters.hh"
#include "G4HepEmTLData.hh"
#include "G4HepEmRandomEngine.hh"

#include "URandom.hh"
#include "Geometry.hh"
#include "PrimaryGenerator.hh"
#include "Results.hh"
#include "EventLoop.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

void SimulationStateDeleter::operator()(const G4HepEmState* state) const noexcept {
  // FreeG4HepEmData releases the nested tables, not their owner.
  if (state->fData) FreeG4HepEmData(state->fData);
  delete state->fData;
  delete state->fParameters;
  delete state;
}

SimulationState load_simulation_state(const std::string& data_file) {
  std::ifstream input(data_file);
  if (!input.is_open())
    throw std::runtime_error("run_simulation: cannot open data file: " + data_file);
  SimulationState state(G4HepEmStateFromJson(input));
  if (!state || !state->fData || !state->fParameters)
    throw std::runtime_error("run_simulation: incomplete state in data file: " + data_file);
  return state;
}

SimulationRunResult simulate(
  const std::string&         data_file,
  const std::vector<double>& absorber_thickness,
  const std::vector<double>& gap_thickness,
  double                     transverse_size,
  const std::string&         particle,
  double                     particle_energy_min,
  double                     particle_energy_max,
  int                        num_events,
  double                     random_seed,
  int                        verbosity,
  int                        gy,
  int                        gz,
  double                     gaussian_sigma,
  const std::vector<double>& bar_edep,
  const std::vector<double>& absorber_thickness_dot,
  const std::vector<double>& gap_thickness_dot,
  bool                       absorber_passive,
  const G4HepEmState*        state)
{
  if (absorber_thickness.size() != gap_thickness.size())
    throw std::runtime_error(
      "run_simulation: absorber_thickness and gap_thickness must have the same length "
      "(got " + std::to_string(absorber_thickness.size()) + " vs "
              + std::to_string(gap_thickness.size()) + ")");

  const int num_layers = (int)absorber_thickness.size();

  if (num_layers == 0)
    throw std::runtime_error("run_simulation: design must contain at least one layer");
  if (num_events < 0)
    throw std::runtime_error("run_simulation: num_events must be non-negative");
  if (gy <= 0 || gz <= 0)
    throw std::runtime_error("run_simulation: gy and gz must be positive");

  if (!absorber_thickness_dot.empty() && (int)absorber_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: absorber_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(absorber_thickness_dot.size()));

  if (!gap_thickness_dot.empty() && (int)gap_thickness_dot.size() != num_layers)
    throw std::runtime_error(
      "run_simulation: gap_thickness_dot must be empty or match the number of layers "
      "(" + std::to_string(num_layers) + "), got " + std::to_string(gap_thickness_dot.size()));

  const std::size_t total_cells_per_layer = static_cast<std::size_t>(gy) * static_cast<std::size_t>(gz);
  const std::size_t total_cells_per_event = static_cast<std::size_t>(num_layers) * total_cells_per_layer;

#ifdef CODI_REVERSE
  const std::size_t num_events_sz = static_cast<std::size_t>(num_events);
  const std::size_t event_cell_size = num_events_sz * total_cells_per_event;
  if (!bar_edep.empty() && bar_edep.size() != total_cells_per_event && bar_edep.size() != event_cell_size) {
    throw std::runtime_error(
      "run_simulation: bar_edep must be empty, per-event (" + std::to_string(total_cells_per_event) +
      "), or per-event-per-cell (" + std::to_string(event_cell_size) + "), got " +
      std::to_string(bar_edep.size()));
  }
#endif

  SimulationRunResult result;

  auto owned_state = state ? SimulationState{} : load_simulation_state(data_file);
  if (!state) state = owned_state.get();

  URandom theURnd(random_seed);
  G4HepEmRandomEngine theRandomEngine(&theURnd);
  G4HepEmTLData theTLData;
  theTLData.SetRandomEngine(&theRandomEngine);

  Geometry theGeometry;
  theGeometry.SetNumLayers(num_layers);
  theGeometry.SetCaloSizeYZ(transverse_size);
  theGeometry.SetGranularity(gy, gz);
  theGeometry.SetGaussianSigma(gaussian_sigma);
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
  if (bar_edep.size() == event_cell_size) {
    theResult.barEdep = bar_edep;
  } else {
    theResult.barEdep.assign(event_cell_size, 0.0);
    if (!bar_edep.empty()) {
      // EventLoop indexes a full run, so repeat a single-event seed for every event.
      for (std::size_t event = 0; event < num_events_sz; ++event)
        std::copy(bar_edep.begin(), bar_edep.end(),
                  theResult.barEdep.begin() + event * total_cells_per_event);
    }
  }
  theResult.barThicknessAbsorber.resize(num_layers);
  theResult.barThicknessGap.resize(num_layers);
#endif

  EventLoop::ProcessEvents(theTLData, *state, thePrimaryGenerator, theGeometry,
                           theResult, num_events, verbosity,
                           particle_energy_min, particle_energy_max, &theURnd);

  result.edep_per_cell = std::move(theResult.fEdepPerCell);
  result.particle_energies = std::move(theResult.fParticleEnergy);
#if CODI_FORWARD
  result.edep_per_cell_dot = std::move(theResult.fEdepPerCellD);
#endif
#ifdef CODI_REVERSE
  result.bar_thickness_absorber.resize(num_layers);
  result.bar_thickness_gap.resize(num_layers);
  for (int i = 0; i < num_layers; ++i) {
    result.bar_thickness_absorber[i] = theResult.barThicknessAbsorber[i].getSum();
    result.bar_thickness_gap[i]      = theResult.barThicknessGap[i].getSum();
  }
  result.bar_particle_energy_mean = theResult.barParticleEnergy.getMean();
  result.bar_particle_energy_var  = theResult.barParticleEnergy.getVar();
#endif

  return result;
}
