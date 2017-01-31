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

class ViewerEditor(AmberBuilder):
    ''' Fixing/Building/Editing/Visualizing in Jupyter notebook

    Inheritance: AmberPDBFixer --> Leapify --> AmberBuilder --> ViewerEditor

    Notes
    -----
    EXPERIMENTAL
    '''

    def __init__(self, *args, **kwargs):
        super(ViewerEditor, self).__init__(*args, **kwargs)
        self._view = None
        # TODO: decorator?
        self.delay_update_structure = False
        self.build_protein = wrap(super(ViewerEditor, self).build_protein, fixer=self)
        self.build_bdna = wrap(super(ViewerEditor, self).build_bdna, fixer=self)
        self.build_adna = wrap(super(ViewerEditor, self).build_adna, fixer=self)
        self.build_arna = wrap(super(ViewerEditor, self).build_arna, fixer=self)
        self.strip = wrap(super(ViewerEditor, self).strip, fixer=self)
        self.add_hydrogen = wrap(super(ViewerEditor, self).add_hydrogen, fixer=self)
        self.add_missing_atoms = wrap(super(ViewerEditor, self).add_missing_atoms, fixer=self)
        self.remove_water = wrap(super(ViewerEditor, self).remove_water, fixer=self)
        self.assign_histidine = wrap(super(ViewerEditor, self).assign_histidine, fixer=self)
        self.pack = wrap(self.pack, fixer=self)
        self.mutate = wrap(self.mutate, fixer=self)
        self.leapify = wrap(self.leapify, fixer=self)

    def visualize(self):
        if self.parm.coordinates.shape[0] == 0:
            self._view = nglview.NGLWidget()
        else:
            self._view = super(ViewerEditor, self).visualize()
        return self._view

    def minimize(self, *args, **kwargs):
        def callback(xyz):
            self._view.coordinates_dict = {0: xyz}
        return super(ViewerEditor, self).minimize(*args, callback=callback, **kwargs)

    def update_structure(self):
        _update_structure(self)
