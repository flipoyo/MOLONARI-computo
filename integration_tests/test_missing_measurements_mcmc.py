"""Integration test of the MCMC with missing temperature measurements (NaN) at a sensor: DREAM and the random walk
compare the model with the measured values only, so the energies stay finite and proposals are accepted. The
column is not solved before the MCMC: the result does not depend on the column's state (DREAM used to need it,
B-33)."""

import warnings

import numpy as np
import pytest

from pyheatmy import Column, Layer
from pyheatmy.config import NSECINDAY, ZERO_CELSIUS
from pyheatmy.synthetic_MOLONARI import synthetic_MOLONARI

from conftest import estimate, one_layer_space

pytestmark = pytest.mark.slow  # every test of this file needs an MCMC run

DEPTH_SENSORS = [0.1, 0.2, 0.3, 0.4]
T_MEAN = ZERO_CELSIUS + 12
SERIES = {
    "offset": 0.0,
    "depth_sensors": DEPTH_SENSORS,
    "param_time_dates": [(2026, 1, 1, 0, 0, 0), (2026, 1, 2, 0, 0, 0), 900],   # 1 day, 96 times
    "param_dH_signal": [[0.1]],
    "param_T_riv_signal": [[T_MEAN], [2.0, NSECINDAY, 0.0]],
    "param_T_aq_signal": [[T_MEAN]],
    "sigma_meas_P": None,
    "sigma_meas_T": None,
    "verbose": False,
}
LAYER = {"name": "sand", "zLow": DEPTH_SENSORS[-1], "IntrinK": 1e-12, "n": 0.1,
         "lambda_s": 2.0, "rhos_cs": 4e6, "q_s": 0.0}
PRIORS = {"IntrinK": ((1e-13, 1e-11), 1e-13), "n": (0.1, 0.0), "lambda_s": (2.0, 0.0),
          "rhos_cs": (4e6, 0.0), "q_s": (0.0, 0.0)}   # one free parameter
N_CELLS = 10
N_ITER = 20
N_BURN_IN_MAX = 3
N_CHAIN_DREAM = 7                                   # at least 2 * delta + 1 for the default delta = 3
ONE_CHAIN = 1                                       # random walk Metropolis
MIDDLE_SENSOR = 1
MISSING_TIMES = [10, 11, 50]


def column_with_gaps(tmp_path):
    series = synthetic_MOLONARI.from_dict(SERIES)
    T_measures = [(date, np.array(t, dtype=float)) for date, t in series._molonariT_data]
    for i in MISSING_TIMES:
        T_measures[i][1][MIDDLE_SENSOR] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        column = Column.from_dict({"river_bed": 1.0, "depth_sensors": DEPTH_SENSORS, "offset": 0.0,
                                   "dH_measures": series._molonariP_data, "T_measures": T_measures,
                                   "nb_cells": N_CELLS, "rac": str(tmp_path)})
    column.set_layers(Layer.from_dict(LAYER))
    return column


@pytest.mark.parametrize("nb_chain", [pytest.param(N_CHAIN_DREAM, id="DREAM"), pytest.param(ONE_CHAIN, id="RWM")])
def test_mcmc_with_missing_measurements(tmp_path, nb_chain):
    column = column_with_gaps(tmp_path)
    assert np.isnan(column.get_measured_temperatures()[MIDDLE_SENSOR, MISSING_TIMES]).all()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # short burn-in
        result = estimate(column, one_layer_space(LAYER["name"], PRIORS), nb_iter=N_ITER, nb_chain=nb_chain,
                          nitmaxburning=N_BURN_IN_MAX)
    assert np.isfinite(result.energies()).all()
    assert np.max(result.acceptance) > 0
