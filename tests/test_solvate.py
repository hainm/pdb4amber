from jamber.builder import solvate, build_protein

def test_solvate():
    traj = build_protein('ALA ALA', ['alpha:1-2'])
    print(traj)
    traj2 = solvate(traj)
    assert traj2.n_atoms == 1092
