from cosmosis.datablock import option_section, names

def setup(options):
    config = {}

    config["reaction_input_section"] = options.get_string(option_section, "reaction_input_section", "reaction")
    config["nonlinear_matter_power_input_section"] = options.get_string(option_section, "nonlinear_matter_power_input_section", names.matter_power_nl)
    config["nonlinear_matter_power_output_section"] = options.get_string(option_section, "nonlinear_matter_power_output_section", names.matter_power_nl)

    return config

def execute(block, config):
    z, k, reaction = block.get_grid(config["reaction_input_section"], "z", "k_h", "reaction")
    z, k, pofk_nonlin = block.get_grid(config["nonlinear_matter_power_input_section"], "z", "k_h", "p_k")

    pofk = reaction*pofk_nonlin
    block.replace_grid(config["nonlinear_matter_power_output_section"], "z", z, "k_h", k, "p_k", pofk)

    return 0

def cleanup(config):
    pass
