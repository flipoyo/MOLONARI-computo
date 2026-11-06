"""Integration test of a prior on IntrinK written in the former -log10(k) notation (11 for 1e-11 m2): the
state space does not know what IntrinK is, but pyheatmy's Layer refuses such values, so the MCMC stops at its first
evaluation of the model, before any state is recorded, with an error that names IntrinK. Moved from the tests of
the priors of pyheatmy's Layer, which held the priors until pyheatmy 0.12.0."""

import pytest

from conftest import DEMO_LAYER, DEMO_PRIORS, estimate, one_layer_space

pytestmark = pytest.mark.slow  # the column of the demo notebook is solved

SIGMA = 0.01
OLD_NOTATION = 11                       # -log10(k) of 1e-11 m2
OLD_NOTATION_RANGE = (4 + 7, 9 + 7)     # prior range of demoPyheatmyMultilayer, in -log10(k)
OLD_NOTATION_PRIORS = [(OLD_NOTATION_RANGE, SIGMA), (OLD_NOTATION, SIGMA)]   # a range and a fixed value
N_CHAIN_DREAM = 7                       # at least 2 * delta + 1 for the default delta = 3
ONE_CHAIN = 1                           # random walk Metropolis
N_ITER = 1
N_BURN_IN_MAX = 1


@pytest.mark.parametrize("prior", OLD_NOTATION_PRIORS)
@pytest.mark.parametrize("nb_chain", [pytest.param(N_CHAIN_DREAM, id="DREAM"), pytest.param(ONE_CHAIN, id="RWM")])
def test_mcmc_stops_on_intrink_in_former_notation(demo_column, prior, nb_chain):
    space = one_layer_space(DEMO_LAYER["name"], {**DEMO_PRIORS, "IntrinK": prior})
    with pytest.raises(ValueError, match="IntrinK"):
        estimate(demo_column, space, nb_iter=N_ITER, nb_chain=nb_chain, nitmaxburning=N_BURN_IN_MAX)
