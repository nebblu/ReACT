import os
import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import pyreact

def setup(options):
    config = {}

    config["module"] = pyreact.ReACT()
    config["verbose"] = options.get_int(option_section, "verbose", 1)
    config["massloop"] = options.get_int(option_section, "massloop", 30)
    config["model"] = options.get_int(option_section, "model", 2)
    config["reaction_output_section"] = options.get_string(option_section, "reaction_output_section", "reaction")
    config["linear_matter_power_output_section"] = options.get_string(option_section, "linear_matter_power_output_section", names.matter_power_lin)

    return config

def execute(block, config):
    h = block[names.cosmological_parameters, "h0"]
    omega_m = block[names.cosmological_parameters, "omega_m"]
    omega_b = block[names.cosmological_parameters, "omega_b"]
    sigma_8 = block[names.cosmological_parameters, "sigma_8"]
    n_s = block[names.cosmological_parameters, "n_s"]
    mg1 = block[names.cosmological_parameters, "mg1"]

    Pk = block[names.matter_power_lin, "p_k"][0]
    k_h = block[names.matter_power_lin, "k_h"]
    z = block[names.matter_power_lin, "z"]

    # np.savetxt("pofk.txt", Pk)
    # np.savetxt("k.txt", k_h)
    # np.savetxt("z.txt", z)

    reaction, pofk_lin = config["module"].compute_reaction(
                                h, n_s, omega_m, omega_b, sigma_8, mg1,
                                z, k_h, Pk, is_transfer=False, mass_loop=config["massloop"], model=config["model"],
                                verbose=config["verbose"])

    block.put_grid(config["reaction_output_section"], "z", z, "k_h", k_h, "reaction", reaction)
    block.replace_grid(config["linear_matter_power_output_section"], "z", z, "k_h", k_h, "p_k", pofk_lin)

    return 0

def cleanup(config):
    pass
