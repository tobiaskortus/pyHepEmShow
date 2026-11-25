#include "ad_type.h"

#ifndef SIMULATIONRUN_HH
#define SIMULATIONRUN_HH

#include <memory>
#include <string>
#include <vector>

struct SimulationRunResult {
  std::vector<double> edep_per_cell;
  std::vector<double> particle_energies;
#if CODI_FORWARD
  std::vector<double> edep_per_cell_dot;
#endif
#ifdef CODI_REVERSE
  std::vector<double> bar_thickness_absorber;
  std::vector<double> bar_thickness_gap;
  double bar_particle_energy_mean = 0.0;
  double bar_particle_energy_var = 0.0;
#endif
};

struct G4HepEmState;
struct SimulationStateDeleter {
  void operator()(const G4HepEmState* state) const noexcept;
};
using SimulationState = std::unique_ptr<const G4HepEmState, SimulationStateDeleter>;
SimulationState load_simulation_state(const std::string& data_file);

SimulationRunResult simulate(
  const std::string& data_file,
  const std::vector<double>& absorber_thickness,
  const std::vector<double>& gap_thickness,
  double transverse_size,
  const std::string& particle,
  double particle_energy_min,
  double particle_energy_max,
  int num_events,
  double random_seed,
  int verbosity,
  int gy,
  int gz,
  double gaussian_sigma,
  const std::vector<double>& bar_edep,
  const std::vector<double>& absorber_thickness_dot,
  const std::vector<double>& gap_thickness_dot,
  bool absorber_passive,
  const G4HepEmState* state = nullptr);

#endif // SIMULATIONRUN_HH
