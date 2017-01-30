from functools import wraps
from .amber_builder import AmberBuilder
import nglview

def _update_structure(fixer):
   if fixer._view is not None:
       traj = nglview.ParmEdTrajectory(fixer.parm)
       if fixer._view.n_components == 0:
           fixer._view.add_trajectory(traj)
           if len(fixer.parm.atoms) <= 50:
               # for small molecule
               fixer._view.add_licorice()
           fixer._view.center()
       else:
           struct = dict(data=traj.get_structure_string(),
                         ext='pdb')
           fixer._view._remote_call('replaceStructure',
                   target='Widget',
                   args=[struct,])

def wrap(func, fixer):
    @wraps(func)
    def me(*args, **kwargs):
        result = func(*args, **kwargs)
        if not fixer.delay_update_structure:
            _update_structure(fixer)
        return result
    return me

class Viewer(AmberBuilder):

    def __init__(self, *args, **kwargs):
        super(Viewer, self).__init__(*args, **kwargs)
        self._view = None
        # TODO: decorator?
        self.delay_update_structure = False
        self.build_protein = wrap(super(Viewer, self).build_protein, fixer=self)
        self.strip = wrap(super(Viewer, self).strip, fixer=self)
        self.add_hydrogen = wrap(super(Viewer, self).add_hydrogen, fixer=self)
        self.add_missing_atoms = wrap(super(Viewer, self).add_missing_atoms, fixer=self)
        self.remove_water = wrap(super(Viewer, self).remove_water, fixer=self)
        self.assign_histidine = wrap(super(Viewer, self).assign_histidine, fixer=self)
        self.pack = wrap(self.pack, fixer=self)
        self.mutate = wrap(self.mutate, fixer=self)
        self.leapify = wrap(self.leapify, fixer=self)

    def visualize(self):
        if self.parm.coordinates.shape[0] == 0:
            self._view = nglview.NGLWidget()
        else:
            self._view = super(Viewer, self).visualize()
        return self._view

    def minimize(self, *args, **kwargs):
        def callback(xyz):
            self._view.coordinates_dict = {0: xyz}
        return super(Viewer, self).minimize(*args, **kwargs, callback=callback)

    def update_structure(self):
        _update_structure(self)
