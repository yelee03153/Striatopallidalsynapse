import numpy as np

# Simulation parameters

dt = 0.001 # Time step for simulation (seconds)
time = np.arange(0, 1, dt) # Time vector from 0 to 1 second with step size dt
spike_times = [0.3, 0.4] # Timing of action potentials (spikes) for paired-pulse stimulation
Ca0 = 50e-9 # Baseline (resting) intracellular calcium concentration (M)
tau_ipsc = 0.02 # Time constant for IPSC exponential decay (seconds)
m = 4 # Hill coefficient for calcium-dependent release probability
n_infty = 1.0 # Steady-state vesicle availability (normalized to 1)

# Observed PPR and conditions
PPR_data = {
    'control': 1,
    '6-OHDA': 0.8,
    'Quinpirole': 2,
    'Quinpirole 6-OHDA': 1.
}

# Dopamine concentrations (M) under each condition
DA_concentrations = {
    'control': 50e-9,
    '6-OHDA': 1e-11,
    'Quinpirole': 50e-9,
    'Quinpirole 6-OHDA': 0.0
}

quinpirole_concentrations = {
    'control': 0.0,
    '6-OHDA': 0.0,
    'Quinpirole': 10e-6,
    'Quinpirole 6-OHDA': 10e-6
}

# Bounds for optimization
param_bounds = [
    (0.1, 0.9),          # P0: Baseline release probability
    (1e-9, 1e-5),        # Kc: Calcium affinity for release (dissociation constant)
    (0.00, 1),           # tau_f: Facilitation decay time constant
    (0.00, 2),           # tau_d: Depression recovery time constant
    (0.00, 2),           # tau_rec: Vesicle replenishment time constant
    (2e-2, 5e-1),        # tau_Ca: Calcium decay time constant
    (1e-7, 100e-5),      # calcium_influx: Total calcium influx per spike
    (0.1, 0.7)           # calcium_channel_open_prob: Opening probability of calcium channels
]
