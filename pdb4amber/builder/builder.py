import warnings
import parmed as pmd

from .compat import StringIO

warnings.filterwarnings("ignore")

def _get_water():
    pdb_str = """
HETATM    1  O   TP3     1       0.000   0.000   0.000  0.00  0.00           O  
HETATM    2  H1  TP3     1       0.957   0.000   0.000  0.00  0.00           H  
HETATM    3  H2  TP3     1      -0.240   0.927   0.000  0.00  0.00           H  
END   
"""
    return pmd.formats.pdb.PDBFile.parse(StringIO(pdb_str))
    
func_dict = {
    'water': _get_water
}

def get(what, *args, **kwargs):
    return func_dict[what](*args, **kwargs)
