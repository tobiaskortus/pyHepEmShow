#include "ad_type.h"


#include "Geometry.hh"

#include "Box.hh"

#include <iostream>
#include <stdexcept>
#include <string>

Geometry::Geometry() {
  // default values: 50 layers of 2.3 [mm] absorber (PbWO4) and 5.7 [mm] gap (lAr)
  fNumLayers  =  50;
  fCaloSizeYZ = 400; // default value [mm]

  // YZ granularity defaults (set explicitly so they are valid before UpdateParameters)
  fCellSizeY    = fCaloSizeYZ;
  fCellSizeZ    = fCaloSizeYZ;
  fInvCellSizeY = 1.0 / fCaloSizeYZ;
  fInvCellSizeZ = 1.0 / fCaloSizeYZ;

  // initialize per-layer thickness vectors with default values
  fAbsThick.assign(fNumLayers, 2.3);
  fGapThick.assign(fNumLayers, 5.7);
  fLayerThick.resize(fNumLayers);
  fLayerStartX.resize(fNumLayers + 1);

  // these will be computed automatically in the `UpdateParameters`
  fCaloStartX       = 0.0;
  fPrimaryXPosition = 0.0;

  // crate shapes here for all objects:
  // - their proper size is set when calling `UpdateParameters` below
  // - material index is set to 0, 1 or 2 that corresponds to (using the default
  //   material list) {0 - G4_Galactic; 1 - G4_PbWO4, 2 - G4_lAr}
  fBoxWorld  = new Box("World", 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxCalo   = new Box("Calo" , 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxLayer  = new Box("Layer", 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxAbs    = new Box("Abs"  , 1, 1,1,1); // material index = 1 G4_PbWO4
  fBoxGap    = new Box("Gap"  , 2, 1,1,1); // material index = 2 G4_lAr

  UpdateParameters();
}


Geometry::~Geometry() {
  delete fBoxWorld;
  delete fBoxCalo;
  delete fBoxLayer;
  delete fBoxGap;
  delete fBoxAbs;
}


void Geometry::UpdateParameters() {
  // compute per-layer thicknesses and cumulative layer start positions
  fLayerThick.resize(fNumLayers);
  fLayerStartX.resize(fNumLayers + 1);
  fLayerStartX[0] = 0.0;
  for (int i = 0; i < fNumLayers; ++i) {
    fLayerThick[i] = fAbsThick[i] + fGapThick[i];
    fLayerStartX[i+1] = fLayerStartX[i] + fLayerThick[i];
  }
  fCaloThick = fLayerStartX[fNumLayers];

  // set/calculate the left hand side x point where the calorimeter starts
  fCaloStartX = -0.5*fCaloThick;
  // set/calculate a world size such that everything fits inside
  const G4double worldThick = 1.1*fCaloThick;
  // set/calculate the mid-point between the world and calorimeter on the left
  fPrimaryXPosition = -0.25*(worldThick + fCaloThick);

  // half size of all (but the world) along the YZ plane
  const double halfCaloYZ = 0.5*fCaloSizeYZ;

  fBoxWorld->SetHalfLength(0.5*worldThick, 0);
  fBoxWorld->SetHalfLength(1.1*halfCaloYZ, 1);
  fBoxWorld->SetHalfLength(1.1*halfCaloYZ, 2);

  int iMatCalo = 0;
  fBoxCalo->SetMaterialIndx(iMatCalo);
  fBoxCalo->SetHalfLength(0.5*fCaloThick, 0);
  fBoxCalo->SetHalfLength(halfCaloYZ, 1);
  fBoxCalo->SetHalfLength(halfCaloYZ, 2);

  // YZ sizes for layer/abs/gap are fixed; x half-lengths are set per-call in CalculateDistanceToOut
  fBoxLayer->SetHalfLength(halfCaloYZ, 1);
  fBoxLayer->SetHalfLength(halfCaloYZ, 2);

  fBoxAbs->SetHalfLength(halfCaloYZ, 1);
  fBoxAbs->SetHalfLength(halfCaloYZ, 2);

  fBoxGap->SetHalfLength(halfCaloYZ, 1);
  fBoxGap->SetHalfLength(halfCaloYZ, 2);

  // set x half-lengths to first-layer values as initial/fallback
  if (fNumLayers > 0) {
    fBoxLayer->SetHalfLength(0.5*fLayerThick[0], 0);
    fBoxAbs->SetHalfLength(0.5*fAbsThick[0], 0);
    fBoxGap->SetHalfLength(0.5*fGapThick[0], 0);
  }

  // update YZ cell sizes from granularity
  fCellSizeY  = fCaloSizeYZ / fGranularityY;
  fCellSizeZ  = fCaloSizeYZ / fGranularityZ;
  fInvCellSizeY = 1.0 / fCellSizeY;
  fInvCellSizeZ = 1.0 / fCellSizeZ;
}


void Geometry::SetAbsThickTangentPerLayer(const std::vector<double>& dots) {
#if CODI_FORWARD
  if (dots.size() != fAbsThick.size())
    throw std::runtime_error(
      "Geometry::SetAbsThickTangentPerLayer: expected " + std::to_string(fAbsThick.size()) +
      " per-layer tangents, got " + std::to_string(dots.size()));
  for (std::size_t i = 0; i < dots.size(); ++i)
    SET_DOTVALUE(fAbsThick[i], dots[i]);
  UpdateParameters();
#else
  (void)dots;
#endif
}


void Geometry::SetGapThickTangentPerLayer(const std::vector<double>& dots) {
#if CODI_FORWARD
  if (dots.size() != fGapThick.size())
    throw std::runtime_error(
      "Geometry::SetGapThickTangentPerLayer: expected " + std::to_string(fGapThick.size()) +
      " per-layer tangents, got " + std::to_string(dots.size()));
  for (std::size_t i = 0; i < dots.size(); ++i)
    SET_DOTVALUE(fGapThick[i], dots[i]);
  UpdateParameters();
#else
  (void)dots;
#endif
}


// note: try to keep this more verbose than fast to keep it clear
G4double Geometry::CalculateDistanceToOut(G4double* r, G4double *v, Box** currentVolume, int* indxLayer, int* indxAbs) {
  // init everything to a step in the `world` case
  *currentVolume = fBoxWorld;
  *indxLayer     = -1;
  *indxAbs       = -1;

  // calculate position in the `calorimeter` system:
  // - only x-coordinate is need as everything is centered along the yz
  // - actually its the same as the global: the calorimeter is not translated nor rotated
  const G4double rx_Calo = r[0];
  const G4double dToCalo = fBoxCalo->DistanceToOut(r, v);
  // check if about leaving the calorimeter volume: distance to out is zero
  if (dToCalo == 0.0) {
    // currentVolume is already set to `world`
    return 1.0E+20;
  }

  // calculate the position in the `layer` system:
  // - find the layer index by scanning cumulative layer start positions
  // - a boundary-limited step can end exactly on an internal boundary, but
  //   floating-point rounding may place the post-step point just before or
  //   just after that boundary. A plain `<` comparison then assigns the
  //   point to the old or the new layer based on that rounding. If it
  //   chooses the old layer, the point lies just outside that layer, so the
  //   stepper takes a push step. That extra step consumes random numbers and
  //   skips the normal physics step, which changes the subsequent shower. To
  //   make the choice independent of rounding, points within `kNavTol` of a
  //   boundary are assigned to the volume the track is entering: the later
  //   layer when moving forward, and the earlier layer when moving backward
  //   or parallel.
  //   This rule is needed for finite-difference comparisons, not because AD
  //   is wrong: forward and reverse AD already agree without it. Finite
  //   differences follow a perturbed trajectory, so without this rule they
  //   can select a different boundary branch and resample the random stream.
  constexpr double kNavTol = 1.0E-9;
  const G4double rxRel = rx_Calo + 0.5*fCaloThick; // position relative to calo left edge
  int iLayer = fNumLayers - 1;
  for (int i = 0; i < fNumLayers - 1; ++i) {
    if (GET_VALUE(rxRel) < GET_VALUE(fLayerStartX[i+1]) - kNavTol ||
        (GET_VALUE(rxRel) <= GET_VALUE(fLayerStartX[i+1]) + kNavTol && GET_VALUE(v[0]) < kNavTol)) {
      iLayer = i;
      break;
    }
  }
  *indxLayer = iLayer;

  // set layer box x half-length for this specific layer
  const G4double layThick = fLayerThick[iLayer];
  const G4double absThick = fAbsThick[iLayer];
  const G4double gapThick = fGapThick[iLayer];
  fBoxLayer->SetHalfLength(0.5*layThick, 0);

  // - then the corresponding translation vector and transform the point to layer-local coords
  const G4double trLayeri = -0.5*fCaloThick + fLayerStartX[iLayer] + 0.5*layThick;
  const G4double rx_Layer = rx_Calo - trLayeri;
  r[0] =  rx_Layer;

  // calculate the distance to the `layer` boundary along the given direction
  // why: tolerance and direction was not considered! So to detect here that
  //      the point is actually miss-located (distance is zero in that case.)
  if (fBoxLayer->DistanceToOut(r, v) == 0.0) {
    return 0.0;
    // NOTE: I could also push here and do recursion but keep it clear and push only in the steppers
  }

  // calculate if the point is in the `absorber` or the `gap` part of the `layer`:
  // - same determinism rule as for the layer scan: on the interface, the
  //   point is assigned to the volume the track is entering
  bool inAbsorber;
  if (gapThick == 0) {
    inAbsorber = true;
  } else {
    const G4double interfaceX = -0.5*layThick + absThick; // interface in layer coords
    if (GET_VALUE(rx_Layer) < GET_VALUE(interfaceX) - kNavTol ||
        (GET_VALUE(rx_Layer) <= GET_VALUE(interfaceX) + kNavTol && GET_VALUE(v[0]) < kNavTol)) {
      inAbsorber = true;
    } else {
      inAbsorber = false;
    }
  }
  if (inAbsorber) { // in the `absorber`
    // calculate the position in the `absorber` system:
    // - the translation vector and transform the point
    fBoxAbs->SetHalfLength(0.5*absThick, 0);
    const G4double trAbs = -0.5*(layThick - absThick);
    r[0] = rx_Layer - trAbs;
    // set what is left and calculate the distance to the `absorber` boundary along
    // the given direction (again, I could push here and do recursion whenever it's zero)
    *currentVolume = fBoxAbs;
    *indxAbs       = 0;
    return fBoxAbs->DistanceToOut(r, v);
  } else { // in the `gap`
    // calculate the position in the `gap` system:
    // - the translation vector and transform the point
    fBoxGap->SetHalfLength(0.5*gapThick, 0);
    const G4double  trGap = -0.5*(layThick - gapThick) + absThick;
    r[0] = rx_Layer - trGap;
    // set what is left and calculate the distance to the `gap` boundary along
    // the given direction (again, I could push here and do recursion whenever it's zero)
    *currentVolume = fBoxGap;
    *indxAbs       = 1;
    return fBoxGap->DistanceToOut(r, v);
  }
}
