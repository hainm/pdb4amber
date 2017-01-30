import os
import parmed as pmd
from pdb4amber.builder.pytraj_build import (build_bdna,
    build_adna, build_arna)
from pdb4amber.builder import nab
from utils import tempfolder

def test_nab():
    fn = 'nuc.pdb'
    with tempfolder():
        nab.run("""
        molecule m; 
        m = fd_helix( "abdna", "aaaaaaaaaa", "dna" );
        putpdb("nuc.pdb", m, "-wwpdb");
        """)

        root = 'tmp_jamber'
        assert not os.path.exists(root + '.nab')
        assert not os.path.exists(root + '.c')
        assert not os.path.exists(root + '.out')

        parm = pmd.load_file(fn)
        assert len(parm.atoms) == 638

def test_nab_bdna():
    with tempfolder():
        parm = build_bdna('AAAAAAAAAA')
        assert len(parm.atoms) == 638

def test_nab_adna():
    with tempfolder():
        parm = build_adna('AAAAAAAAAA')
        assert len(parm.atoms) == 638

def test_nab_arna():
    with tempfolder():
        parm = build_arna('AAAAAAAAAA')
        print(parm.atoms)
