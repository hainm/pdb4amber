from pdb4amber import builder

def test_water():
    water = builder.get('water')
    assert len(water.residues) == 1
    assert water.residues[0].name == 'TP3'
