from .pdb4amber import AmberPDBFixer
from .utils import tempfolder
from .leap_runner import run_tleap


class Leapify(AmberPDBFixer):

    def minimize(self, **kwargs):
        ''' port ParmED'minimize function

        Examples
        --------
        >>> # Using PME
        >>> fixer.minimize(igb=None, saltcon=None, cutoff=10., maxcyc=100, tol=1E-6)
        '''
        igb = kwargs.get('igb', None)

        if igb is not None:
            # GB
            if 'saltcon' not in kwargs:
                kwargs['saltcon'] = 0.
            if 'cutoff' not in kwargs:
                kwargs['cutoff'] = 999.
        else:
            # PME
            kwargs['saltcon'] = None
            if 'cutoff' not in kwargs:
                kwargs['cutoff'] = 8.

        if 'maxcyc' not in kwargs:
            kwargs['maxcyc'] = 500.
        if 'tol' not in kwargs:
            kwargs['tol'] = 1E-6

        from parmed.tools.simulations import sanderapi
        sanderapi.minimize(self.parm, **kwargs)

    def leapify(self, *args, **kwargs):
        with tempfolder():
            ns_names = self.find_non_starndard_resnames() if self.parm is not None else []
            gaplist = self.find_gaps() if self.parm is not None else []
            sslist, _ = self.find_disulfide() if self.parm is not None else []
            self.parm = run_tleap(self.parm, ns_names, gaplist, sslist,
                    *args, **kwargs)
