from pdb4amber.builder.pytraj_build import solvate, build_protein

def test_solvate():
    parm = build_protein('ALA ALA', ['alpha:1-2'])
    print(parm)
    parm2 = solvate(parm)
    assert len(parm2.atoms) == 1092
