import sys
sys.path.append("../../pyreact")

import numpy as np
import react

def test_reaction_module():
    mod = react.ReactionModule()

    h = 0.68
    n_s = 0.9645
    Omega_m = 0.308489
    Omega_b = 0.0475779
    sigma_8 = 0.815088
    fR0 = 1e-5
    mass_loop = 30

    # k = np.loadtxt("k.txt")
    # z = np.loadtxt("z.txt")
    # pofk = np.loadtxt("pofk.txt")
    #
    # reaction, p_lin = mod.compute_reaction(h, n_s, Omega_m, Omega_b, sigma_8, fR0, mass_loop,
    #                                 z, k, pofk, is_transfer=False)

    k, t = np.loadtxt("../benchmarks/transfer.dat", unpack=True)
    z = np.loadtxt("../benchmarks/reaction_z.dat")

    reaction, p_lin = mod.compute_reaction(h, n_s, Omega_m, Omega_b, sigma_8, fR0, mass_loop,
                                    z, k, t, is_transfer=True)

    reaction_target = np.loadtxt("../benchmarks/reaction.dat")

    assert np.allclose(reaction, reaction_target)

if __name__ == "__main__":
    test_reaction_module()
