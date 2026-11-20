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

## Native HDF5 CLI

A standalone C++ executable can read a design and write flat HDF5 simulation
outputs for DSBO's native file contract. It uses the header-only `dsbo::io` HDF5
SDK and requires HDF5 >= 2.0 with thread safety enabled.

Build the forward (non-AD) and reverse-mode variants with the helper script:

```sh
bash externals/configure_native.sh forward /path/to/dsbo/cpp/include /path/to/hdf5
bash externals/configure_native.sh backward /path/to/dsbo/cpp/include /path/to/hdf5
```

Each build produces `hepemshow/build_native_<mode>/HepEmShowNative`. The script
initializes the G4HepEm install for the selected mode, disables the Python
extension, and enables the native target. It accepts an optional build type as
the fourth argument (`Release` by default).

The executable takes `forward` or `backward` as its first argument. It reads the
design from an HDF5 TensorDict with `absorber` and `gap` fields:

```sh
HepEmShowNative forward \
  --input in.h5 \
  --edep-output simulation.h5 \
  --energy-output energies.h5 \
  --labels-output labels.h5 \
  --data-file /opt/hepemshow/data/hepem_data_fe_pbwo4.json \
  --num-events 1000 --seed 1234 --particle e- \
  --energy-min 10000 --energy-max 10000 \
  --transverse-size 400 --gy 1 --gz 1 \
  --absorber-passive --verbosity 0
```

Forward mode writes three flat primal files: `simulation.h5` with shape
`(num_events, num_layers, gy, gz)`, and `energies.h5` and `labels.h5` with shape
`(num_events,)`. The executable does not write the design vector; callers such
as the DSBO container module write it from the design tensors they already hold.

Reverse mode reads the flat edep primal and an incoming adjoint, reruns the
simulation with reverse-mode CoDiPack using the same settings and seed, and
writes a TensorDict gradient with `absorber` and `gap` keys:

```sh
HepEmShowNative backward \
  --input in.h5 \
  --edep-output simulation.h5 \
  --edep-adjoint simulation_gradient.h5 \
  --input-adjoint in_gradient.h5 \
  --data-file /opt/hepemshow/data/hepem_data_fe_pbwo4.json \
  --num-events 1000 --seed 1234 --particle e- \
  --energy-min 10000 --energy-max 10000 \
  --transverse-size 400 --gy 1 --gz 1 \
  --absorber-passive --verbosity 0
```

The backward command is only available when G4HepEm was built with
`CODI_REVERSE`; otherwise it reports an error and exits nonzero. Forward-mode
tangents are not supported by the native CLI.

## Usage

```python
import pyhepemshow

settings = dict(
    data_file=pyhepemshow.get_hepemshow_data_file(),
    absorber_thickness=[2.3] * 5,
    gap_thickness=[5.7] * 5,
    particle="e-",
    particle_energy_min=10_000.0,
    particle_energy_max=10_000.0,
    num_events=10,
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
- `gy` and `gz` set transverse cell counts (default: `1` each).
- `absorber_passive=True` excludes absorber energy deposition from the readout.

All modes return a dictionary with these NumPy arrays:

| Key | Shape | Meaning |
| --- | --- | --- |
| `edep_per_cell` | `(num_events, num_layers, gy, gz)` | Energy deposited in each cell |
| `particle_energies` | `(num_events,)` | Incident particle kinetic energies |

Materials and physics tables are defined by `data_file`. Different materials
require a new file generated with
[HepEmShow-DataGeneration](https://hepemshow.readthedocs.io/en/latest/IntroAndInstall/components.html#the-hepemshow-datageneration-application-main),
which requires Geant4.

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
`edep_per_cell`.

### Reverse Mode

`reverse` computes a vector-Jacobian product. Supply the objective's gradient
with respect to energy deposition as `bar_edep`. For total deposited energy:

```python
result_ad = pyhepemshow.reverse.run_simulation(
    **settings,
    bar_edep=[1.0] * result["edep_per_cell"].size,
)
absorber_gradient = result_ad["bar_thickness_absorber"]
gap_gradient = result_ad["bar_thickness_gap"]
```

Supply `bar_edep` for all events and cells, flattened in C order. Both returned
thickness gradients have shape `(num_layers,)`.

## Reference

[1] M. Aehle, M. Novák, V. Vassilev, N. R. Gauger, L. Heinrich, M. Kagan, D. Lange.
[Optimization using pathwise algorithmic derivatives of electromagnetic shower simulations](https://doi.org/10.1016/j.cpc.2024.109491).
*Computer Physics Communications* **309**, 109491 (2025).

## License and Attribution

The original [HepEmShow](https://github.com/SciCompKL/hepemshow) is released under
the Apache License 2.0. CoDiPack is released under GPLv3; conveying the combined
work requires compliance with the GPL.
