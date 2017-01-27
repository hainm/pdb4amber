from functools import wraps
from .leapify import Leapify
from nglview import ParmEdTrajectory

def _update_structure(fixer):
   if fixer._view is not None:
       traj = ParmEdTrajectory(fixer.parm)
       struct = dict(data=traj.get_structure_string(),
                     ext='pdb')
       fixer._view._remote_call('replaceStructure',
               target='Widget',
               args=[struct,])

def wrap(func, fixer):
    @wraps(func)
    def me(*args, **kwargs):
        result = func(*args, **kwargs)
        _update_structure(fixer)
        return result
    return me

class Viewer(Leapify):

    def __init__(self, *args, **kwargs):
        super(Viewer, self).__init__(*args, **kwargs)
        self._view = None
        self.strip = wrap(super(Viewer, self).strip, fixer=self)
        self.add_hydrogen = wrap(super(Viewer, self).add_hydrogen, fixer=self)
        self.add_missing_atoms = wrap(super(Viewer, self).add_missing_atoms, fixer=self)
        self.remove_water = wrap(super(Viewer, self).remove_water, fixer=self)
        self.assign_histidine = wrap(super(Viewer, self).assign_histidine, fixer=self)
        self.pack = wrap(self.pack, fixer=self)
        self.mutate = wrap(self.mutate, fixer=self)
        self.leapify = wrap(self.leapify, fixer=self)

    def visualize(self):
        self._view = super(Viewer, self).visualize()
        return self._view

    def minimize(self, *args, **kwargs):
        def callback(xyz):
            self._view.coordinates_dict = {0: xyz}
        return super(Viewer, self).minimize(*args, **kwargs, callback=callback)
