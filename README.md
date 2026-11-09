# MOLONARI-computo

The computing part of MOLONARI: estimating the exchanges of water and heat between a river and its bed from the
temperatures and the head difference measured in the river bed. This repository is the root of the MOLONARI-computo
tree, in which cgitsync places the following repositories:

| Path | Repository | Content |
|---|---|---|
| `sources/pyheatmy` | [MOLONARI-computo-pyheatmy](https://github.com/flipoyo/MOLONARI-computo-pyheatmy) | `pyheatmy`, the direct model: heat and water transport in a 1D column of the river bed (`Column`, `Layer`, `Column.simulate`) |
| `sources/bayesian_inversion` | [MOLONARI-computo-bayesian-inversion](https://github.com/flipoyo/MOLONARI-computo-bayesian-inversion) | `bayesian_inversion`, the MCMC (DREAM, random walk Metropolis) that estimates the parameters of a model from its measurements; for a column, `run_mcmc(column.simulate, ...)` |
| `tuto` | [MOLONARI-computo-tuto](https://github.com/flipoyo/MOLONARI-computo-tuto) | tutorials: a synthetic data set generated with `pyheatmy`, its Bayesian inversion, and the inversion of the real data of `Point034` |

The two packages do not import each other; the tests that use them together are in `integration_tests/` here.
Molonaviz 2 (MOLONARI-API) uses both of them as references for the tests of its own engine.

## Installation

The tree is described by `install_molonari_computo.cgs`, at the root of this repository. With
[ComplexGitSync](https://github.com/flipoyo/ComplexGitSync), run from its own clone:

```sh
cgitsync bootstrap <path to>/install_molonari_computo.cgs
```

clones the four repositories at their place in a new workspace, on their branch `CGS-MOLO`. Then, with
[pixi](https://pixi.sh), at the root of the tree:

```sh
pixi install
pixi run test-all
```

`pixi.toml` installs the two packages from `sources/` in editable mode. `pixi run test-all` runs the tests of
`pyheatmy`, then those of `bayesian_inversion`, then those of `integration_tests/`. The tutorials
run with the Python of `.pixi/envs/default` as the kernel (see `tuto/README.md`). Each package can also be installed
on its own with pip; its tests run on GitHub Actions at every push to its repository.

## This repository

- `integration_tests/`: the tests of the two packages together, with copies of the data they read, run with
  `pixi run test-shared`.
- `data/`: processed data of real monitoring points (`Point021`, `Point034`).
- `research/`: notebooks of the MOLONARI 2023 to 2025 editions, kept as a record of their work and not maintained;
  they were written for former versions of `pyheatmy`. `adimensionnement/`: the equations in non-dimensional form;
  `multilayer/`: a column of several layers; `physics/`: the dependence of the parameters on temperature, the
  propagation of the error of the pressure measurement and the remanence time.
- `TD4_analyse_donnees.ipynb`: a practical session on the processing of real sensor data before the direct model;
  it needs raw data that are not in this repository, and is to move to MOLONARI-data.

The former `pyheatmy` package, version 0.9.7, the demonstration notebooks written for it, the descriptions of the
sensors and shafts (`configuration/`), the LoRa emitter sample and the other research projects were in this
repository until the release `MOLONARI-computo-2` of the tree; the tutorials of `tuto/` replace the notebooks, and
all of them are in the history of this repository, for instance at the tag `MOLONARI-computo-1`.

## License

See [LICENSE](LICENSE).

2021 Mathis Bourdin & Youri Tchouboukoff
2022 Amélie Impéror, Antoine Poirier, Guillaume de Rochefort, Loris Megy, Paul Bonin Ciosi, Valentin Alleaume, Xinbei Jiang
2023 Mattéo Leturcq-Daligaux, Nicolas Matte, Mathis Chevé, Zhan Jing, Dan Maurel

Project manager : Nicolas Flipo
