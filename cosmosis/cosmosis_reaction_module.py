import os
import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import pyreact

def setup(options):
    config = {}

    config["module"] = pyreact.ReACT()
    config["mode"] = options.get_string(option_section, "mode", "fR").lower()
    if config["mode"] != "fr":
        raise ValueError(f"ReACT mode {config['mode']} not supported.")

    config["verbose"] = options.get_int(option_section, "verbose", 1)
    config["massloop"] = options.get_int(option_section, "massloop", 30)
    config["z_max"] = options.get_double(option_section, "z_max", 2.5)
    config["reaction_output_section"] = options.get_string(option_section, "reaction_output_section", "reaction")
    config["linear_matter_power_output_section"] = options.get_string(option_section, "linear_matter_power_output_section", names.matter_power_lin)

    return config

def execute(block, config):
    h = block[names.cosmological_parameters, "h0"]
    omega_m = block[names.cosmological_parameters, "omega_m"]
    omega_b = block[names.cosmological_parameters, "omega_b"]
    sigma_8 = block[names.cosmological_parameters, "sigma_8"]
    n_s = block[names.cosmological_parameters, "n_s"]
    if config["mode"] == "fr":
        mg1 = block[names.cosmological_parameters, "fR0"]

    Pk = block[names.matter_power_lin, "p_k"]
    k_h = block[names.matter_power_lin, "k_h"]
    z = block[names.matter_power_lin, "z"]

    z_react = z[z < config["z_max"]]

    reaction, pofk_lin = config["module"].compute_reaction(
                                h, n_s, omega_m, omega_b, sigma_8, mg1, 
                                z_react, k_h, Pk[0], is_transfer=False, mass_loop=config["massloop"],
                                verbose=config["verbose"])

    # Replace linear power spectrum below z_max with MG linear power spectrum
    Pk[z < config["z_max"]] = pofk_lin
    # Pad the reaction with 1 above z_max
    reaction = np.concatenate((reaction, np.ones((np.count_nonzero(z >= config["z_max"]), len(k_h)))), axis=0)

    block.put_grid(config["reaction_output_section"], "z", z, "k_h", k_h, "reaction", reaction)
    block.replace_grid(config["linear_matter_power_output_section"], "z", z, "k_h", k_h, "p_k", Pk)

    return 0

def cleanup(config):
    pass
