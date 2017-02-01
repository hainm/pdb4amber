# adapted from NGLVIEW
# https://github.com/arose/nglview

import sys
import unittest
import warnings
warnings.filterwarnings('ignore')
PY3 = sys.version_info[0] == 3

from pdb4amber.visualizer import ViewerEditor

from utils import get_fn

try:
    import scipy
    has_scipy = True
except ImportError:
    has_scipy = False

try:
    from ipywidgets import Widget
    from ipykernel.comm import Comm
    import nglview
    has_nglview = True

    #------------------------------------------------------
    # Utility stuff from ipywidgets tests: create DummyComm
    # we dont need Jupyter notebook for testing
    #------------------------------------------------------
    class DummyComm(Comm):
        comm_id = 'a-b-c-d'
    
        def open(self, *args, **kwargs):
            pass
    
        def send(self, *args, **kwargs):
            pass
    
        def close(self, *args, **kwargs):
            pass
    
    _widget_attrs = {}
    displayed = []
    undefined = object()

    _widget_attrs['_comm_default'] = getattr(Widget, '_comm_default', undefined)
    Widget._comm_default = lambda self: DummyComm()
    _widget_attrs['_ipython_display_'] = Widget._ipython_display_
    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()
    Widget._ipython_display_ = raise_not_implemented
except ImportError:
    has_nglview = False
    nglview = Comm = DummyComm = _widget_attrs = displayed = undefined = Widget  = None

@unittest.skipUnless(has_nglview and has_scipy, 'require nglview, ipywidgets and scipy')
def test_viewer_editor():
    """
    """
    editor = ViewerEditor(get_fn('2igd/2igd.pdb'))
    view = editor.visualize()
    view
    editor.remove_water()
    editor.leapify()
    editor.minimize(igb=None, maxcyc=10, disp=False)
