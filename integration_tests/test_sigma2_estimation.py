"""Integration tests of the estimation of the temperature measurement variance (bayesian_inversion.run_mcmc with
typealgo="estimesigma", on pyheatmy's Column.simulate). The layer parameters are fixed at the values that generated the synthetic data, so
sigma2 is the only unknown and its posterior is known in closed form: with n measurements, a sum of squared
residuals S, the Gaussian likelihood (sigma2)^(-n/2) exp(-S / (2 sigma2)) and the Jeffreys prior 1 / sigma2,
sigma2 follows an inverse gamma distribution of shape n/2 and scale S/2, of mean S / (n - 2)."""

import numpy as np
import pytest

from pyheatmy import Column, Layer
from pyheatmy.config import NSECINDAY, ZERO_CELSIUS
from pyheatmy.synthetic_MOLONARI import synthetic_MOLONARI

from bayesian_inversion import mcmc
from bayesian_inversion.config import SIGMA2_MAX_T, SIGMA2_MIN_T

from conftest import SEED, estimate, fixed_priors, one_layer_space

pytestmark = pytest.mark.slow  # every test of this file needs an MCMC run

DEPTH_SENSORS = [0.1, 0.2, 0.3, 0.4]
T_MEAN = ZERO_CELSIUS + 12
SIGMA_T = 0.1                                     # K, so sigma2 = 0.01, inside the prior range (0.001, 1)
SERIES = {
    "offset": 0.0,
    "depth_sensors": DEPTH_SENSORS,
    "param_time_dates": [(2026, 1, 1, 0, 0, 0), (2026, 1, 11, 0, 0, 0), 900],   # 10 days, 960 times
    "param_dH_signal": [[0.1]],
    "param_T_riv_signal": [[T_MEAN], [2.0, NSECINDAY, 0.0]],
    "param_T_aq_signal": [[T_MEAN]],
    "sigma_meas_P": None,
    "sigma_meas_T": SIGMA_T,
    "verbose": False,
}
LAYER = {"name": "sand", "zLow": DEPTH_SENSORS[-1], "IntrinK": 1e-12, "n": 0.1,
         "lambda_s": 2.0, "rhos_cs": 4e6, "q_s": 0.0}
FIXED_SPACE = one_layer_space(LAYER["name"], fixed_priors(LAYER))
N_CELLS = 20

N_BURN_IN_MAX = 10
N_CHAIN_DREAM = 7                                 # at least 2 * delta + 1 for the default delta = 3
N_ITER_DREAM = 60
ONE_CHAIN = 1                                     # random walk Metropolis
N_ITER_RWM = 150
# 2880 measurements: the posterior standard deviation of sigma2 is sqrt(2 / (n - 4)) ~ 2.6 % of its mean
MEAN_RTOL = 0.03

N_ITER_CONTRACT = 20
N_BURN_IN_CONTRACT = 5


def column_with_fixed_layer(dH_measures, T_measures, tmp_path):
    column = Column.from_dict({"river_bed": 1.0, "depth_sensors": DEPTH_SENSORS, "offset": 0.0,
                               "dH_measures": dH_measures, "T_measures": T_measures,
                               "nb_cells": N_CELLS, "rac": str(tmp_path)})
    column.set_layers(Layer.from_dict(LAYER))
    return column


@pytest.fixture(scope="module")
def measurements(tmp_path_factory):
    """Noisy sensor temperatures generated with the layer parameters, as in demo_genData."""
    np.random.seed(SEED)
    series = synthetic_MOLONARI.from_dict(SERIES)
    generating = Column.from_dict({"river_bed": 1.0, "depth_sensors": DEPTH_SENSORS, "offset": 0.0,
                                   "dH_measures": series._molonariP_data, "T_measures": series._T_Shaft_measures,
                                   "nb_cells": N_CELLS, "rac": str(tmp_path_factory.mktemp("generating"))})
    generating.set_layers(Layer.from_dict(LAYER))
    series._measures_column_one_layer(generating)
    return series._molonariP_data, series._molonariT_data


@pytest.fixture(scope="module")
def inverse_gamma_mean(measurements, tmp_path_factory):
    """Posterior mean S / (n - 2) of sigma2, from the residuals of the direct model with the fixed parameters."""
    column = column_with_fixed_layer(*measurements, tmp_path_factory.mktemp("direct"))
    simulated = column.simulate([{}])["sensors"]
    measured = column.get_measured_temperatures()
    n = measured.size
    return np.sum((measured - simulated) ** 2) / (n - 2)


def run_estimesigma(measurements, tmp_path, **mcmc_arguments):
    column = column_with_fixed_layer(*measurements, tmp_path)
    return estimate(column, FIXED_SPACE, typealgo="estimesigma", **mcmc_arguments)


@pytest.mark.parametrize("nb_chain, nb_iter", [pytest.param(N_CHAIN_DREAM, N_ITER_DREAM, id="DREAM"),
                                             pytest.param(ONE_CHAIN, N_ITER_RWM, id="RWM")])
def test_sigma2_posterior_mean_matches_inverse_gamma(measurements, inverse_gamma_mean, nb_chain, nb_iter, tmp_path):
    result = run_estimesigma(measurements, tmp_path, nb_iter=nb_iter, nb_chain=nb_chain, nitmaxburning=N_BURN_IN_MAX)
    sigma2 = result.sigma2()
    second_half = sigma2[len(sigma2) // 2:]       # the chains start anywhere in the prior range
    assert np.mean(second_half) == pytest.approx(inverse_gamma_mean, rel=MEAN_RTOL)


@pytest.mark.parametrize("nb_chain", [pytest.param(N_CHAIN_DREAM, id="DREAM"), pytest.param(ONE_CHAIN, id="RWM")])
def test_energy_uses_physical_sigma2_and_flat_prior_in_mcmc_space(measurements, nb_chain, tmp_path, monkeypatch):
    """Every energy is computed with a variance inside the prior range, and the prior term -log(density) is the
    same for all of them: sigma2 is uniform in log10(sigma2), which is the Jeffreys prior 1 / sigma2."""
    calls = []
    original_energy = mcmc.compute_energy

    def recording_energy(temp_simul, temp_ref, sigma2, sigma2_distrib):
        calls.append((float(sigma2), float(sigma2_distrib(sigma2))))
        return original_energy(temp_simul, temp_ref, sigma2, sigma2_distrib)

    monkeypatch.setattr(mcmc, "compute_energy", recording_energy)
    result = run_estimesigma(measurements, tmp_path, nb_iter=N_ITER_CONTRACT, nb_chain=nb_chain,
                             nitmaxburning=N_BURN_IN_CONTRACT)
    sigma2, density = np.array(calls).T
    assert np.all((SIGMA2_MIN_T <= sigma2) & (sigma2 <= SIGMA2_MAX_T))
    assert np.all(density == density[0])
    recorded = result.sigma2()
    assert np.all((SIGMA2_MIN_T <= recorded) & (recorded <= SIGMA2_MAX_T))
