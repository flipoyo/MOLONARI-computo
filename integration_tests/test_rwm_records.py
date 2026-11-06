"""Integration tests of what the random walk Metropolis records (bayesian_inversion.run_mcmc with one chain) on the
demo column: the states and the space they are stored in, the physical parameters of the result, the acceptance
rate and the quantiles."""

import numpy as np
import pytest

pytestmark = pytest.mark.slow  # every test of this file needs an MCMC run

# Settings of the reference RWM run
N_ITER = 100
N_BURN_IN_MAX = 20
ONE_CHAIN = 1

ENERGY_RTOL = 1e-6
QUANTILE_RTOL = 1e-6


@pytest.fixture(scope="module")
def rwm_result(run_demo_mcmc):
    """One RWM run, shared by the tests of this file."""
    return run_demo_mcmc(nb_iter=N_ITER, nb_chain=ONE_CHAIN, nitmaxburning=N_BURN_IN_MAX)


def state_after(states, iteration):
    """State after a main loop iteration: the history starts with the best initial state."""
    return states[iteration + 1]


def test_rwm_records_initial_state_and_one_state_per_iteration(rwm_result):
    assert len(rwm_result.states) == N_ITER + 1


def test_every_rwm_state_energy_matches_its_parameters_in_mcmc_space(rwm_result, replay):
    """States are stored in the MCMC space, as DREAM stores them: converting them to physical values and solving
    the direct model gives back the recorded energy."""
    unique_states = {id(state): state for state in rwm_result.states}.values()  # rejections repeat a state
    inconsistent = [state.energy for state in unique_states
                    if replay(state)[0] != pytest.approx(state.energy, rel=ENERGY_RTOL)]
    assert not inconsistent, f"{len(inconsistent)} of {len(unique_states)} states, energies {inconsistent[:5]}"


def test_rwm_physical_parameters_are_in_prior_range(rwm_result, demo_space):
    """best and params convert the states to physical values inside the prior range."""
    intrinK = demo_space.param_names.index("IntrinK")
    lower, upper = demo_space.priors[0][intrinK].user_range
    best = rwm_result.best()[0]["IntrinK"]
    every = rwm_result.params()[:, 0, intrinK]
    assert lower <= best <= upper
    assert np.all((lower <= every) & (every <= upper))


def test_rwm_acceptance_rate_is_accepted_over_proposed(rwm_result):
    """The rate is the fraction of iterations where the state changed."""
    states = rwm_result.states
    changes = [not np.array_equal(state_after(states, i - 1).params, state_after(states, i).params)
               for i in range(N_ITER)]
    assert rwm_result.acceptance == pytest.approx(np.mean(changes))


def test_rwm_quantiles_use_only_accepted_states(rwm_result, replay):
    """The quantiles are those of the fields of the accepted states.
    With these settings every iteration is kept for the quantiles."""
    fields = [replay(state_after(rwm_result.states, i))[1:] for i in range(N_ITER)]
    temperatures = np.array([t for t, _ in fields], np.float32)
    flows = np.array([f for _, f in fields], np.float32)

    for quantile, values in rwm_result.quantiles["temperatures"].items():
        np.testing.assert_allclose(values, np.quantile(temperatures, quantile, axis=0), rtol=QUANTILE_RTOL)
    for quantile, values in rwm_result.quantiles["flows"].items():
        np.testing.assert_allclose(values, np.quantile(flows, quantile, axis=0), rtol=QUANTILE_RTOL)
