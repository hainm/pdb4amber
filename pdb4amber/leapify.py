import subprocess
from .pdb4amber import AmberPDBFixer
from parmed.tools.simulations import sanderapi
from .utils import tempfolder

class Leapify(AmberPDBFixer):

    def minimize(self, *args, **kwargs):
        sanderapi.minimize(self.parm, *args, **kwargs)

    def leapify(self):
        with tempfolder():
            r
