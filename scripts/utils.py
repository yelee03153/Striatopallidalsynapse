# Functions for the synptic modeling 
#
#
#
#
import numpy as np

base_DA = 60e-9 #50- 100 nM Dodson P. D., et al. & Magill P. J. Representation of spontaneous movement by dopaminergic neurons is cell-type selective and disrupted in parkinsonism. PNAS (2016)
# D2 autoreceptors have been shown to inhibit the P/Q and N-type channels (Cardazo & Bean, 1995).
def D2R_activation(DA_concentration, agonist_concentration=0, EC50_DA=10e-9, EC50_agonist=20e-9, max_activation=0.5):
    """
    Calculates the D2 receptor activation based on dopamine and agonist concentrations.
    """
    # Activation due to dopamine
    activation_DA =1* DA_concentration / (EC50_DA + DA_concentration)
    
    # Activation due to quinpirole (agonist)
    activation_agonist = agonist_concentration / (EC50_agonist + agonist_concentration)
    
    # Total activation (ensure it does not exceed max_activation)
    total_activation = max_activation * (activation_DA + activation_agonist)
    total_activation = np.clip(total_activation, 0, max_activation)
    
    return total_activation

def is_spike_time(t, spike_times, dt):
    return any(np.isclose(t, spike_times, atol=dt/2))

