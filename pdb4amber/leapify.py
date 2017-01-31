from .pdb4amber import AmberPDBFixer
from .utils import tempfolder, which
from .leap_runner import run_tleap
import subprocess


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

    def _run_antechamber_ccif(self, residue_name, use_am1_and_maxcyc_zero=True):
        # Original code from amber_adaptbx
        '''
        run antechamber from a components.cif file:
        '''
        if not which('antechamber'):
            raise OSError("Make sure to have antechamber")
        mask = ':' + residue_name
        ccif = 'test.mol2'
        self.parm[mask][':1'].save(ccif, overwrite=True)
        cmds = []
        cmd = 'antechamber -i %s -fi mol2 -bk %s -o %s.mol2 -fo mol2 \
          -s 2 -pf y -c bcc -at gaff2' % (ccif, residue_name, residue_name)
        if use_am1_and_maxcyc_zero:
            cmd += ' -ek "qm_theory=\'AM1\', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"'
        cmds.append(cmd)

        for cmd in cmds:
            output = subprocess.check_output(cmd, shell=True)

        cmd = 'parmchk2 -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' % (residue_name,
                                                                 residue_name)
        subprocess.check_call(cmd, shell=True)
