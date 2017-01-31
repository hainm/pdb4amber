from pdb4amber.amber_builder import AmberBuilder

def test_solvate():
    builder = AmberBuilder()
    builder.build_protein('ALA ALA', ['alpha:1-2'])
    builder.solvate()
    assert len(builder.parm.atoms) == 1092
