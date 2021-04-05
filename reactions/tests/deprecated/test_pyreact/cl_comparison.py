import numpy as np

import cosmosis.runtime.config
import cosmosis.runtime.pipeline
import cosmosis.datablock

import pyccl as ccl

pi = np.pi

PIPELINE_FILE = "pipeline.ini"

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "omega_k" : 0.0,
                                               "w"       : -1.0,
                                               "mnu"     : 0.06,
                                               "fR0"     : 1.0e-5},
                  "halo_model_parameters"   : {"A"       : 2.5,
                                               "eta0"    : 0.603,}}


    pipeline_ini = cosmosis.runtime.config.Inifile(PIPELINE_FILE)
    values_ini = cosmosis.runtime.config.Inifile(None)
    values_ini.read_dict(param_dict)

    reaction_pipeline = cosmosis.runtime.pipeline.LikelihoodPipeline(pipeline_ini, values=values_ini)
    reaction_data = reaction_pipeline.run_parameters([])
    
    modules = pipeline_ini.get("pipeline", "modules").split(" ")
    modules.remove("reaction")
    modules.remove("multiply_reaction")
    pipeline_ini.set("pipeline", "modules", " ".join(modules))

    fiducial_pipeline = cosmosis.runtime.pipeline.LikelihoodPipeline(pipeline_ini, values=values_ini)
    fiducial_data = fiducial_pipeline.run_parameters([])

    # Plot results

    k = reaction_data["matter_power_nl", "k_h"]
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(5,5))
    fig.subplots_adjust(hspace=0, left=0.15)

    ax[0].loglog(k, reaction_data["matter_power_nl", "p_k"][0], label="Reaction")
    ax[0].loglog(k, fiducial_data["matter_power_nl", "p_k"][0], label="No reaction")
    ax[1].semilogx(k, reaction_data["matter_power_nl", "p_k"][0]/fiducial_data["matter_power_nl", "p_k"][0] - 1)

    ax[0].legend(frameon=False)
    ax[0].set_ylabel("$P(k)$ [$h^{-3}$ Mpc$^{3}$]")
    # ax[1].set_ylim(-0.05, 0.05)
    ax[1].set_xlabel("$k$ [$h$ Mpc$^{-1}$]")
    ax[1].set_ylabel("$\Delta P(k)/P(k)$ [$h$ Mpc$^{-1}$]")
    
    fig.suptitle("Reaction, P(k)")
    fig.savefig("reaction_pofk.pdf")

    n_bin = reaction_data["shear_cl", "nbin"]
    ell = reaction_data["shear_cl", "ell"]
    u = ell**2

    fig, ax = plt.subplots(n_bin, n_bin, sharex=True, sharey=True, figsize=(7, 5))
    fig.subplots_adjust(hspace=0, wspace=0)

    ell_range = 10, 10000

    for i in range(n_bin):
        for j in range(n_bin):
            if j > i:
                ax[i][j].axis("off")
            else:
                b = f"bin_{i+1}_{j+1}"
                ax[i][j].loglog(ell, u*reaction_data["shear_cl", b], label="Reaction")
                ax[i][j].loglog(ell, u*fiducial_data["shear_cl", b], label="No reaction")

                ax[i][j].set_xlim(*ell_range)
                
    ax[0,0].legend(fontsize="small", frameon=False, loc=2, bbox_to_anchor=(1,1))

    for p in ax[-1]:
        p.set_xlabel(r"$\ell$")
    for p in ax:
        p[0].set_ylabel(r"$\ell^2\ C_\ell$")

    fig.suptitle("Reaction, Cls")
    fig.savefig("reaction_cl.pdf")
    plt.show()