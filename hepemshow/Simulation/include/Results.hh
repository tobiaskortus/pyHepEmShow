#include "ad_type.h"


#ifndef RESULTS_HH
#define RESULTS_HH

/**
 * @file    Results.hh
 * @author  M. Novak
 * @date    July 2023
 *
 * @brief A collection of data that are recorded during the simulation.
 *
 * The following data is recorded during the simulation (mean is per event):
 *  - mean values in the individual layers of the calorimeter for energy deposit,
 *    neutral (gamma) and charged (electron/positron) particle simulation steps
 *  - mean number of energy deposited in the `absorber` and `gap`
 *  - mean number of secondary gamma, electron and positrons produced
 *  - mean number of neutral (gamma) and charged (electron/positron)
 *
 * Quantities, recorded in the individual layers are stored in histograms and
 * written to files at the end of the simulation while the others are reported
 * in the screen. An example looks like
 * ```
 *       --- Results::WriteResults ----------------------------------
 *
 *       Absorber: mean Edep = 6722.95 [MeV] and  Std-dev = 309.636 [MeV]
 *       Gap     : mean Edep = 2571.75 [MeV] and  Std-dev = 118.507 [MeV]
 *
 *       Mean number of gamma       4457.043
 *       Mean number of e-          7957.899
 *       Mean number of e+          428.922
 *
 *       Mean number of e-/e+ steps 36097
 *       Mean number of gamma steps 40436.2
 *       ------------------------------------------------------------
 * ```
 */

#include "Hist.hh"
#include <vector>
#include "accumulator.hh"

/**
 * Data that are collected during the entire `run` of the simulation.
 * Mean quantities are computed over the simulated events.
 */
struct Results {
  Hist fEdepPerLayer_CurrentEvent;  ///< per-layer energy deposit accumulator for the current event
  std::vector<double> fEdepPerLayer; ///< per-event per-layer energy deposits, flat row-major [events*layers]
  std::vector<double> fParticleEnergy; ///< per-event primary particle kinetic energy [MeV]
  #ifdef CODI_FORWARD
    std::vector<double> fEdepPerLayerD; ///< per-event per-layer dot values, flat row-major [events*layers]
  #endif
  #ifdef CODI_REVERSE
    std::vector<double> barEdep; ///< Bar values of the edeps, to be set in the beginning.
    std::vector<Accumulator<double>> barThicknessAbsorber, barThicknessGap; ///< Per-layer accumulated bar values of the absorber/gap thicknesses.
    Accumulator<double> barParticleEnergy; ///< Accumulated bar value of the particle energy.
    std::vector<G4double> pThicknessAbsorber, pThicknessGap; ///< Per-layer copies of the thickness variables used as AD inputs.
    G4double pParticleEnergy; ///< Copy of the particle energy variable used as AD input.
  #endif
};

#endif // RESULTS_HH
