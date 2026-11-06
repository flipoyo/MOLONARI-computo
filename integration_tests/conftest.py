"""Shared fixtures of the tests of the MOLONARI-computo tree, which use pyheatmy and bayesian_inversion together:
the Point034 data set and the homogeneous column of the demo notebook, its state space, an MCMC run on it, and the
replay of a recorded state through the direct model. pytest loads this file automatically."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyheatmy import Column, Layer
from pyheatmy.config import ZERO_CELSIUS

from bayesian_inversion import Prior, StateSpace, compute_energy, run_mcmc

TEST_DATA = Path(__file__).resolve().parent / "data"  # copies of data/ and configuration/ of MOLONARI-computo
DATA_DIR = TEST_DATA / "Point034_processed"
PRESSURE_SENSOR_FILE = TEST_DATA / "P508.csv"
PRESSURE_SENSOR = "P508"
DATE_FORMAT = "%y/%m/%d %H:%M"
TEMPERATURE_COLUMNS = ["temperature_10", "temperature_20", "temperature_30", "temperature_40"]

SEED = 0  # seed of the reference runs of the refactoring
PARAMETERS = ("IntrinK", "n", "lambda_s", "rhos_cs", "q_s")  # pyheatmy's Param

# Demo notebook: column, homogeneous layer and priors
DEMO_COLUMN = {
    "river_bed": 1.0,
    "depth_sensors": [0.1, 0.2, 0.3, 0.4],
    "offset": 0.0,
    "sigma_meas_P": None,
    "sigma_meas_T": None,
    "inter_mode": "lagrange",
    "nb_cells": 40,
}
DEMO_LAYER = {"name": "Layer_homogenous", "zLow": 0.4, "IntrinK": 1e-12, "n": 0.1,
              "lambda_s": 2.0, "rhos_cs": 4e6, "q_s": 0.0}
DEMO_PRIORS = {
    "IntrinK": ((1e-18, 1e-15), 5e-17),
    "n": ((0.01, 0.25), 0.0125),
    "lambda_s": ((1, 10), 0.5),
    "rhos_cs": ((1e6, 1e7), 9e5),
    "q_s": ((-1e-6, 1e-6), 1e-11),
}


def uniform_density(_sigma2):
    """Flat prior on sigma2, as run_mcmc uses in the "no sigma" mode."""
    return 1.0


def one_layer_space(layer_name, priors):
    """State space of a column of one layer: priors maps each parameter of Param to the arguments of its Prior."""
    return StateSpace({layer_name: {name: Prior(*priors[name]) for name in PARAMETERS}})


def fixed_priors(layer):
    """Every parameter fixed at the value of the layer."""
    return {name: (layer[name], 0.0) for name in PARAMETERS}


def estimate(column, space, **mcmc_arguments):
    """run_mcmc on the direct model of the column, from the same seed. Prior uses the random module and DREAM uses
    numpy: run_mcmc seeds both."""
    return run_mcmc(column.simulate, column.get_measured_temperatures(), space, verbose=False, seed=SEED,
                    **mcmc_arguments)


@pytest.fixture(scope="session")
def point034_measures():
    """River head and temperature, and sensor temperatures, processed as in the demo notebook."""
    river = pd.read_csv(DATA_DIR / "processed_pressures.csv",
                        names=["dates", "tension", "temperature_river"], skiprows=1)
    sensors = pd.read_csv(DATA_DIR / "processed_temperatures.csv",
                          names=["dates", *TEMPERATURE_COLUMNS], skiprows=1)
    calibration = pd.read_csv(PRESSURE_SENSOR_FILE, index_col=0)[PRESSURE_SENSOR]
    intercept = float(calibration["Intercept"])
    du_dh = float(calibration["dU/dH"])
    du_dt = float(calibration["dU/dT"])

    river["dates"] = pd.to_datetime(river["dates"], format=DATE_FORMAT)
    sensors["dates"] = pd.to_datetime(sensors["dates"], format=DATE_FORMAT)
    head = (river["tension"] - intercept - river["temperature_river"] * du_dt) / du_dh
    river_temperature = river["temperature_river"] + ZERO_CELSIUS
    sensor_temperatures = sensors[TEMPERATURE_COLUMNS].to_numpy() + ZERO_CELSIUS

    dH_measures = list(zip(river["dates"], zip(head, river_temperature)))
    T_measures = list(zip(sensors["dates"], sensor_temperatures))
    return dH_measures, T_measures


@pytest.fixture(scope="session")
def demo_column(point034_measures, tmp_path_factory):
    """The demo column with its layer. Output files of the column go to a temporary folder, not to the home
    directory. simulate leaves the column unchanged, so the tests share it."""
    dH_measures, T_measures = point034_measures
    column = Column.from_dict({
        **DEMO_COLUMN,
        "dH_measures": dH_measures,
        "T_measures": T_measures,
        "rac": str(tmp_path_factory.mktemp("column_output")),
    })
    column.set_layers(Layer.from_dict(DEMO_LAYER))
    return column


@pytest.fixture(scope="session")
def demo_space():
    return one_layer_space(DEMO_LAYER["name"], DEMO_PRIORS)


@pytest.fixture(scope="session")
def run_demo_mcmc(demo_column, demo_space):
    """Factory: each call runs the MCMC with the given arguments on the demo column, from the same seed."""
    def run(**mcmc_arguments):
        return estimate(demo_column, demo_space, **mcmc_arguments)

    return run


@pytest.fixture(scope="session")
def replay(demo_column, demo_space):
    """Function: solves the direct model with the parameters of a state, and returns its energy, temperatures and
    flows."""
    def solve(state):
        output = demo_column.simulate(demo_space.to_physical(state.params))
        energy = compute_energy(output["sensors"], demo_column.get_measured_temperatures(), state.sigma2,
                                uniform_density)
        return energy, np.asarray(output["temperatures"]), np.asarray(output["flows"])

    return solve


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: runs an MCMC, takes about a minute (skip with -m 'not slow')"
    )
