import numpy as np

def derive_s8_fixed_alpha(chain):
    alpha=0.50
    #alpha=0.77
    omega_m=chain['cosmological_parameters--omega_m']
    try:
        sigma_8=chain['cosmological_parameters--sigma8_input']
    except:
        try:
            sigma_8=chain['cosmological_parameters--sigma_8']
        except:
            try:
                sigma_8=chain['COSMOLOGICAL_PARAMETERS--SIGMA_8']
            except: #sigma8 isn't in the chain; is entirely unconstrained
                # return array of random values in the range [0.5,1.5]
                nentries = omega_m.size
                sigma_8 = np.random.rand(nentries) + 0.5
                
    s8 = sigma_8 * (omega_m/0.3)**alpha
    return s8,"cosmological_parameters--s8_%.3f"%alpha