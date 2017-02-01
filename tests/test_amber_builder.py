import unittest
import parmed as pmd
from numpy.testing import assert_almost_equal as aa_eq
from pdb4amber.amber_builder import AmberBuilder
from pdb4amber.utils import tempfolder

from utils import get_fn

try:
    import pytraj
except ImportError:
    pytraj = None

def _compute_major_groove(builder):
    with  tempfolder():
        pdb_out = 'out.pdb'
        builder.write_pdb(pdb_out)
        traj = pytraj.load(pdb_out)

        data = pytraj.nastruct(traj)
        return data.major[1][0]


@unittest.skipUnless(pytraj is not None, "require pytraj")
def test_solvate():
    builder = AmberBuilder()
    builder.build_protein('ALA ALA', ['alpha:1-2'])
    builder.solvate()
    assert len(builder.parm.atoms) == 1092

@unittest.skipUnless(pytraj is not None, "require pytraj")
def test_build_nuleic_acid():
    builder = AmberBuilder()
    builder.build_bdna('GGGGGG')
    aa_eq(_compute_major_groove(builder), [ 17.24558067])
    aa_eq(_compute_major_groove(builder.build_adna('GGGGGG')), [ 15.72079277])

    # build_xxx will repalce old structure
    aa_eq(_compute_major_groove(builder.build_arna('GGGGGG')), [ 15.18226528])
    with tempfolder():
        pdb_out = 'out.pdb'
        builder.write_pdb(pdb_out)
        with open(pdb_out) as fh:
            assert "O2'" in fh.read()

def test_build_unitcell():
    pdb_fn = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fn)
    builder = AmberBuilder(pdb_fn)
    builder.build_unitcell()
    assert parm.coordinates.shape[0] * 4 == builder.parm.coordinates.shape[0]
