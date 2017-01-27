from .pdb4amber import AmberPDBFixer
from parmed.tools.simulations import sanderapi
from .utils import tempfolder
from .leap_runner import run_tleap


class Leapify(AmberPDBFixer):
    def minimize(self, *args, **kwargs):
        sanderapi.minimize(self.parm, *args, **kwargs)

    def leapify(self):
        with tempfolder():
            ns_names = self.find_non_starndard_resnames()
            gaplist = self.find_gaps()
            sslist, _ = self.find_disulfide()
            self.parm = run_tleap(self.parm, ns_names, gaplist, sslist)
