# pyHepEmShow

pyHepEmShow provides Python bindings for
[HepEmShow](https://hepemshow.readthedocs.io/en/latest/), a compact electromagnetic
shower simulator for layered calorimeters. It supports per-layer geometry
configuration and forward- and reverse-mode automatic differentiation for
detector optimization, with simulation results returned as NumPy arrays.

The bindings use [nanobind](https://github.com/wjakob/nanobind), with
[G4HepEm](https://github.com/SciCompKL/g4hepem) for particle transport and
[CoDiPack](https://github.com/tobiaskortus/CoDiPack) for automatic differentiation.

## Installation

Requires Python 3.10+, a C++17 compiler, CMake, Make, Bash, and Git. On macOS,
install the Xcode command-line tools; on Linux, `patchelf` is needed when
bundling shared libraries.

```sh
git clone --recurse-submodules https://github.com/tobiaskortus/pyHepEmShow.git
cd pyHepEmShow
python -m pip install numpy
python -m pip install -e .
```

The build compiles the simulation with and without automatic differentiation,
using the bundled G4HepEm and CoDiPack submodules.
Geant4 is not required to run simulations.
For an existing checkout, initialize dependencies with
`git submodule update --init --recursive` before installing.

## Usage

```python
import pyhepemshow

settings = dict(
    data_file=pyhepemshow.get_hepemshow_data_file(),
    absorber_thickness=[2.3] * 5,
    gap_thickness=[5.7] * 5,
    particle="e-",
    particle_energies=[10_000.0] * 10,
    random_seed=1234,
    verbosity=0,
)

result = pyhepemshow.forward.run_simulation(**settings)
print(result["edep_per_cell"].shape)  # (10, 5, 1, 1)
print(result["particle_energies"].shape)  # (10,)
```

### Inputs and Outputs

- Thicknesses and `transverse_size` are in mm; energies are in MeV.
- Absorber and gap thickness lists require one entry per layer and equal lengths.
- `particle` accepts `"e-"`, `"e+"`, or `"gamma"`.
- `particle_energies` gives the primary kinetic energy for each event; its length
  determines `num_events`.
- `gy` and `gz` set transverse cell counts (default: `1` each).
- `kernel` selects the transverse readout kernel: `"hard"`, `"gaussian"`, or
  `"trilinear"` (default: `"trilinear"`).
- `gaussian_sigma` sets the Gaussian smearing width; negative values use half a
  cell size.
- `absorber_passive=True` excludes absorber energy deposition from the readout.

All modes return a dictionary with these NumPy arrays:

| Key | Shape | Meaning |
| --- | --- | --- |
| `edep_per_cell` | `(num_events, num_layers, gy, gz)` | Energy deposited in each cell |
| `particle_energies` | `(num_events,)` | Incident particle kinetic energies |

Materials and physics tables are defined by `data_file`. The bundled
`hepem_data.json` uses PbWO4 for the absorber (material index 1) and liquid
argon for the gap (material index 2). Different materials require a new file
generated with
[HepEmShow-DataGeneration](https://hepemshow.readthedocs.io/en/latest/IntroAndInstall/components.html#the-hepemshow-datageneration-application-main),
which requires Geant4; its current default material list is Fe and PbWO4.

## Automatic Differentiation

The AD modes compute pathwise derivatives with respect to layer thicknesses
[[1]](#reference). The examples below reuse `settings` and `result` from above.
Use identical settings and random seeds for simulation and differentiation.

### Forward Mode

`forward_ad` computes a Jacobian-vector product. To differentiate each output
with respect to the first absorber thickness:

```python
absorber_tangent = [0.0] * len(settings["absorber_thickness"])
absorber_tangent[0] = 1.0

result_ad = pyhepemshow.forward_ad.run_simulation(
    **settings,
    absorber_thickness_dot=absorber_tangent,
)
edep_derivative = result_ad["edep_per_cell_dot"]
```

`absorber_thickness_dot` and `gap_thickness_dot` contain one tangent per layer;
omitted tangents are zero. `edep_per_cell_dot` has the same shape as
`edep_per_cell`. For per-event tangents, pass `absorber_thickness_dots` or
`gap_thickness_dots` as `num_events x num_layers` arrays, and
`particle_energies_dot` as a per-event vector. Each event then uses its own
seed stream, so derivatives can be combined across events.

### Reverse Mode

`reverse` computes a vector-Jacobian product. Supply the objective's gradient
with respect to energy deposition as `bar_edep`. For total deposited energy:

```python
import numpy as np

result_ad = pyhepemshow.reverse.run_simulation(
    **settings,
    bar_edep=np.ones_like(result["edep_per_cell"]),
)
absorber_gradient = result_ad["bar_thickness_absorber"]
gap_gradient = result_ad["bar_thickness_gap"]
particle_energy_gradient = result_ad["bar_particle_energies"]
```

`bar_edep` must have shape `(num_events, num_layers, gy, gz)`. The returned
thickness gradients have shape `(num_layers,)`, and `bar_particle_energies` has
shape `(num_events,)`.

## Reference

[1] M. Aehle, M. Novák, V. Vassilev, N. R. Gauger, L. Heinrich, M. Kagan, D. Lange.
[Optimization using pathwise algorithmic derivatives of electromagnetic shower simulations](https://doi.org/10.1016/j.cpc.2024.109491).
*Computer Physics Communications* **309**, 109491 (2025).

## License and Attribution

The original [HepEmShow](https://github.com/SciCompKL/hepemshow) is released under
the Apache License 2.0. CoDiPack is released under GPLv3; conveying the combined
work requires compliance with the GPL.

## Technical notes

- [Navigation tolerance at layer boundaries](docs/navigation_tolerance.md):
  why the boundary tie-break is needed for stable AD gradients and finite
  differences.
