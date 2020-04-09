import os

import numpy as np
import pyreact

HERE = os.path.abspath(os.path.dirname(__file__))

def test_reaction_module():
    mod = pyreact.ReACT()

    h = 0.68
    n_s = 0.9645
    Omega_m = 0.308489
    Omega_b = 0.0475779
    sigma_8 = 0.815088
    fR0 = 1e-5
    mass_loop = 30

    k = np.loadtxt(os.path.join(HERE, "k.txt"))
    z = np.loadtxt(os.path.join(HERE, "z.txt"))
    pofk = np.loadtxt(os.path.join(HERE, "pofk.txt"))

    # Check that the code runs with power spectra
    reaction, p_lin = mod.compute_reaction(h, n_s, Omega_m, Omega_b, sigma_8, fR0,
                                    z, k, pofk, is_transfer=False, mass_loop=mass_loop, verbose=1)

    k, t = np.loadtxt(os.path.join(HERE, "../benchmarks/transfer.dat"), unpack=True)
    z = np.loadtxt(os.path.join(HERE, "../benchmarks/reaction_z.dat"))

    # Check with transfer function and compare with benchmark.
    reaction, p_lin = mod.compute_reaction(h, n_s, Omega_m, Omega_b, sigma_8, fR0,
                                    z, k, t, is_transfer=True, mass_loop=mass_loop, verbose=1)
    
    reaction_target = np.loadtxt(os.path.join(HERE, "../benchmarks/reaction.dat"))
    assert np.allclose(reaction, reaction_target)

if __name__ == "__main__":
    test_reaction_module()
