import parmed
from .leapify import Leapify
from .utils import tempfolder, which, easy_call

class AmberBuilder(Leapify):
    ''' Require many programs in AmberTools (pytraj, tleap, nab, ...)

    Inheritance: AmberPDBFixer --> Leapify --> AmberBuilder
    '''

    def build_protein(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import build_protein
        self.parm = build_protein(*args, **kwargs)
        return self

    def build_bdna(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import build_bdna
        self.parm = build_bdna(*args, **kwargs)
        return self

    def build_adna(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import build_adna
        self.parm = build_bdna(*args, **kwargs)
        return self

    def build_arna(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import build_arna
        self.parm = build_bdna(*args, **kwargs)
        return self

    def solvate(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import solvate
        self.parm = solvate(self.parm, *args, **kwargs)
        return self

    def build_unitcell(self):
        '''

        Requires
        --------
        UnitCell program (AmberTools)
        '''
        UnitCell = which('UnitCell')
        if self.parm.box is None or self.parm.symmetry is None:
            raise ValueError("Must have symmetry and box data")
            raise OSError("Can not find UnitCell program")

        with tempfolder():
            inp_pdb = 'inp.pdb'
            out_pdb = 'out.pdb'
            self.parm.save(inp_pdb)
            out = easy_call([
                'UnitCell',
                '-p', inp_pdb,
                '-o', out_pdb,
            ])
            self.parm = parmed.load_file(out_pdb)
