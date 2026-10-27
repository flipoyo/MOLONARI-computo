import numpy as np
from pyheatmy.config import *

from pyheatmy.params import Param, Prior

class Layer:
    def __init__(
        self,
        name: str,
        zLow: float,  # profondeur en m du bas de la couche (>0)
        IntrinK: float = INTRINK_DEFAULT,
        n: float = N_DEFAULT,
        lambda_s: float = LAMBDA_S_DEFAULT,
        rhos_cs: float = RHOS_CS_DEFAULT,
        q_s: float = Q_S_DEFAULT,
        Prior_IntrinK: Prior = None, # Prior sur IntrinK, c'est un tuple ex: (range=(9, 15), sigma= 0.1)
        Prior_n: Prior = None,
        Prior_lambda_s: Prior = None,
        Prior_rhos_cs: Prior = None,
        Prior_q_s: Prior = None,
    ):
        # Voir si tout ça ne peut pas être simplifié avec un @dataclass

        self.name = name
        self.zLow = zLow
       
        # On stocke les priors
        self.Prior_IntrinK = Prior_IntrinK
        self.Prior_n = Prior_n
        self.Prior_lambda_s = Prior_lambda_s
        self.Prior_rhos_cs = Prior_rhos_cs
        self.Prior_q_s = Prior_q_s
        self.Prior_list = [Prior_IntrinK, Prior_n, Prior_lambda_s, Prior_rhos_cs, Prior_q_s]

        # On crée un tuple de paramètres physiques initiaux
        physical_params = Param(IntrinK, n, lambda_s, rhos_cs, q_s)
        
        # On crée l'attribut mcmc_params en traduisant les valeurs physiques initiales
        # vers l'espace de travail de la MCMC..
        mcmc_vals = [
            # Si un prior existe, on l'utilise pour traduire la valeur.
            prior.physical_to_mcmc(val) if prior is not None else val
            # Sinon (si prior est None), on considère que la valeur MCMC est la même que la valeur physique.
            for prior, val in zip(self.Prior_list, physical_params)
        ]
        self.mcmc_params = Param(*mcmc_vals)

    def sample(self):
        """Met à jour mcmc_params avec des valeurs initiales tirées des priors."""
        sampled_vals = [prior.sample() for prior in self.Prior_list]
        self.mcmc_params = Param(*sampled_vals)
    
    def perturb(self):
        """Met à jour mcmc_params en appliquant une perturbation à chaque valeur."""
        perturbed_vals = [
            prior.perturb(val)
            for prior, val in zip(self.Prior_list, self.mcmc_params)
        ]
        self.mcmc_params = Param(*perturbed_vals)
        
    def __repr__(self) -> str:
        return self.name + f" : ends at {self.zLow} m. " + self.mcmc_params.__repr__()
    
    def set_priors_from_dict(self, priors_dict):
        """Assigne les priors depuis un dictionnaire."""
        # Les valeurs physiques se lisent avec les priors actuels, avant de les remplacer.
        physical_params = self.get_physical_params()
        self.Prior_IntrinK = Prior(*priors_dict['Prior_IntrinK'])
        self.Prior_n = Prior(*priors_dict['Prior_n'])
        self.Prior_lambda_s = Prior(*priors_dict['Prior_lambda_s'])
        self.Prior_rhos_cs = Prior(*priors_dict['Prior_rhos_cs'])
        self.Prior_q_s = Prior(*priors_dict['Prior_q_s'])
        self.Prior_list = [self.Prior_IntrinK, self.Prior_n, self.Prior_lambda_s, self.Prior_rhos_cs, self.Prior_q_s]
        
        # Il faut aussi mettre à jour les mcmc_params lorsque les priors sont assignés !
        mcmc_vals = [
            prior.physical_to_mcmc(val) if prior is not None else val
            for prior, val in zip(self.Prior_list, physical_params)
        ]
        self.mcmc_params = Param(*mcmc_vals)
    
    @classmethod
    def from_dict(cls, monolayer_dict):
        return cls(**monolayer_dict)

    def get_physical_params(self, mcmc_params=None):
        """
        Retourne les paramètres en valeurs PHYSIQUES.
        Cette version est maintenant tolérante à l'absence de Priors.
        """
        if mcmc_params is None:
            mcmc_params = self.mcmc_params

        # CORRECTION : On rend la traduction inverse conditionnelle.
        physical_vals = [
            # Si un prior existe, on l'utilise pour traduire la valeur.
            prior.mcmc_to_physical(val) if prior is not None else val
            # Sinon (si prior est None), on considère que la valeur est déjà physique.
            for prior, val in zip(self.Prior_list, mcmc_params)
        ]
        return Param(*physical_vals)


def getListParameters(layersList, z_cells):
    """
    Renvoie les profils des paramètres physiques sur les cellules de la colonne.
    z_cells : profondeurs des centres des cellules (Column._z_solve), mesurées depuis le lit de la rivière.
    Une cellule appartient à la couche telle que zLow_prev < z <= zLow ; la dernière couche
    couvre toutes les cellules restantes (prolongée ou coupée au fond de la colonne).
    """
    listParameters = np.zeros((len(z_cells), len(Param._fields)))
    zLow_prev = 0
    for i, layer in enumerate(layersList):
        if i == len(layersList) - 1:
            mask = z_cells > zLow_prev
        else:
            mask = (z_cells > zLow_prev) & (z_cells <= layer.zLow)
        listParameters[mask, :] = layer.get_physical_params()
        zLow_prev = layer.zLow
    return tuple(listParameters[:, i] for i in range(len(Param._fields)))


# Add-on function to create the list of layers.
def layersListCreator(layers_spec):
    """
    Construct and return a list of Layer instances.

    Accepts:
    - layers_spec: iterable of tuples or dicts describing each layer.
      If an element is a dict, Layer.from_dict is used.
      If an element is a sequence, it is passed to Layer(*sequence).

    Expected tuple order (compatible with existing calls):
      (name, zLow, IntrinK, n, lambda_s, rhos_cs[, q_s])

    Returns:
      list of Layer
    """
    layers = []
    for spec in layers_spec:
        if isinstance(spec, dict):
            layer = Layer.from_dict(spec)
        else:
            # allow lists/tuples with the exact arguments for Layer.__init__
            layer = Layer(*spec)
        layers.append(layer)
    return layers
