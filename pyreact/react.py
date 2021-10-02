import os
import ctypes as ct
import numpy as np

def array_ctype(ndim, dtype=np.float64, flags="C"):
    return [ct.POINTER(ct.c_int)]*ndim + [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)]

def array_arg(a):
    arr = a
    return (*(ct.c_int(s) for s in arr.shape), arr)

class ReACT:
    libname = "libreact_wrapper.so"
    module_name = "reaction_module"

    def __init__(self):
        self.load_lib()

    def load_lib(self, path=None):
        if path is None:
            path = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(path, self.libname))
        self.lib = ct.CDLL(libpath)

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def test_func(self, a):
        f = self.get_function("test_func")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=2, dtype=np.float64)]

        r = f(*array_arg(a))
        return r

# Compute reaction without massive neutrinos
    def compute_reaction(self, h, n_s, omega_m, omega_b, sigma_8,
                               z, k, Pk,
                               model="f(r)", fR0=None, Omega_rc=None, w=None, wa=None,
                               is_transfer=False, mass_loop=30,
                               verbose=True):


        if max(z) > 2.5:
            raise ValueError("ReACT is unstable above z=2.5, try limiting the range of z values.")

        if len(z) > 1 and z[0] > z[-1]:
            raise ValueError("The z array needs to be ordered from low to high redshift.")

        if model.lower() == "f(r)":
            reaction_model = 2
            modified_gravity_param = fR0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "dgp":
            reaction_model = 3
            modified_gravity_param = Omega_rc
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "gr":
            reaction_model = 1
            modified_gravity_param = 0.0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "quintessence":
            reaction_model = 4
            modified_gravity_param = w
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "cpl":
            reaction_model = 5
            modified_gravity_param = w
            modified_gravity_param2 = wa
            modified_gravity_param3 = 0.0
        else:
            raise ValueError(f"model '{model}' not supported.")

        if modified_gravity_param is None:
            raise ValueError("fR0, Omega_rc or w0 need to be specified.")

        f = self.get_function("compute_reaction")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # P(k, z=0)
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # sigma_8
                      ct.POINTER(ct.c_double),     # model parameter 1 : for f(R) this is fr0, for dgp this is Omega_rc, for CPL or quintessence this is w0
                      ct.POINTER(ct.c_double),     # model parameter 2 : for CPL this is wa
                      ct.POINTER(ct.c_double),     # model parameter 3
                      ct.POINTER(ct.c_int),        # mass_loop
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                      np.ctypeslib.ndpointer(ndim=1, dtype=np.float64, flags="C"),     # modified sigma_8 storage variable
                      ct.POINTER(ct.c_int),        # verbose
                     ]
        reaction = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        sigma8 = np.zeros(1, dtype=Pk.dtype, order="C")

        r =   f(*array_arg(np.ascontiguousarray(Pk, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)), # ReACT expect z order from high to low
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(sigma_8),
                ct.c_double(modified_gravity_param),
                ct.c_double(modified_gravity_param2),
                ct.c_double(modified_gravity_param3),
                ct.c_int(mass_loop),
                ct.c_int(reaction_model),
                *array_arg(reaction),
                *array_arg(p_lin),
                sigma8,
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")
        # Get into CAMB ordering (z, k), with increasing z
        #Tilman: to check output of modsig8
        return reaction[:,::-1].T, p_lin[:,::-1].T, sigma8[0]




    def compute_reaction_nu(self, h, n_s, omega_m, omega_b, omega_nu, As,
                               z, k, Tm, Tcb, klcdm, Tcblcdm,
                               pscale = 0.05,
                               model="f(r)", fR0=None, Omega_rc=None, w=None, wa=None,
                               is_transfer=False, mass_loop=30,
                               verbose=True):

        if max(z) > 2.5:
            raise ValueError("ReACT is unstable above z=2.5, try limiting the range of z values.")

        if len(z) > 1 and z[0] > z[-1]:
            raise ValueError("The z array needs to be ordered from low to high redshift.")

        if model.lower() == "f(r)":
            reaction_model = 2
            modified_gravity_param = fR0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "dgp":
            reaction_model = 3
            modified_gravity_param = Omega_rc
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "gr":
            reaction_model = 1
            modified_gravity_param = 0.0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "quintessence":
            reaction_model = 4
            modified_gravity_param = w
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "cpl":
            reaction_model = 5
            modified_gravity_param = w
            modified_gravity_param2 = wa
            modified_gravity_param3 = 0.0
        else:
            raise ValueError(f"model '{model}' not supported.")

        if modified_gravity_param is None:
            raise ValueError("fR0, Omega_rc or w0 need to be specified.")


        f = self.get_function("compute_reaction_nu")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # full matter transfer
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      *array_ctype(ndim=1, dtype=np.float64), # CDM + baryon transfer
                      *array_ctype(ndim=1, dtype=np.float64), # LCDM CDM+baryon transfer
                      *array_ctype(ndim=1, dtype=np.float64), # k for LCDM transfer
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # omega_nu
                      ct.POINTER(ct.c_double),     # A_s
                      ct.POINTER(ct.c_double),     # pscale
                      ct.POINTER(ct.c_double),     # model parameter 1 : for f(R) this is fr0, for dgp this is Omega_rc, for CPL or quintessence this is w0
                      ct.POINTER(ct.c_double),     # model parameter 2 : for CPL this is wa
                      ct.POINTER(ct.c_double),     # model parameter 3 : for DS this is xi * h
                      ct.POINTER(ct.c_int),        # mass_loop
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                      np.ctypeslib.ndpointer(ndim=1, dtype=np.float64, flags="C"),     # modified sigma_8 storage variable
                      ct.POINTER(ct.c_int),        # verbose
                     ]

        reaction = np.zeros((len(k), len(z)), dtype=k.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=k.dtype, order="C")
        sigma8 = np.zeros(1, dtype=k.dtype, order="C")


        r =   f(*array_arg(np.ascontiguousarray(Tm, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)), # ReACT expect z order from high to low
                *array_arg(np.ascontiguousarray(Tcb, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(Tcblcdm, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(klcdm, dtype=np.float64)),
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(omega_nu),
                ct.c_double(As), ct.c_double(pscale),
                ct.c_double(modified_gravity_param),
                ct.c_double(modified_gravity_param2),
                ct.c_double(modified_gravity_param3),
                ct.c_int(mass_loop),
                ct.c_int(reaction_model),
                *array_arg(reaction),
                *array_arg(p_lin),
                sigma8,
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")

        return reaction[:,::-1].T, p_lin[:,::-1].T, sigma8[0]
