from .leapify import Leapify
from nglview import ParmEdTrajectory

def wrap(func, fixer):
    def me(*args, **kwargs):
        result = func(*args, **kwargs)
        if fixer._view is not None:
            traj = ParmEdTrajectory(fixer.parm)
            struct = dict(data=traj.get_structure_string(),
                          ext='pdb')
            fixer._view._remote_call('replaceStructure',
                    target='Widget',
                    args=[struct,])
        return result
    return me

class Viewer(Leapify):

    def __init__(self, *args, **kwargs):
        super(Viewer, self).__init__(*args, **kwargs)
        self._view = None
        self.strip = wrap(super(Viewer, self).strip, fixer=self)
        self.add_hydrogen = wrap(super(Viewer, self).add_hydrogen, fixer=self)

    def visualize(self):
        self._view = super(Viewer, self).visualize()
        return self._view

