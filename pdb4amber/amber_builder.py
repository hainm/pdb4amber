from .leapify import Leapify

class AmberBuilder(Leapify):
    def build_protein(self, *args, **kwargs):
        from pdb4amber.builder.pytraj_build import build_protein
        self.parm = build_protein(*args, **kwargs)
        return self
