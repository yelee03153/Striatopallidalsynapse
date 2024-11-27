import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Define the model parameters and time constants
Ca0 = 30e-9
Kd_DA = 10e-7  # Dissociation constant for D2 receptor (10 nM converted to M)
n_Hill = 2  # Hill coefficient for DA effect on calcium
P0 = 0.5
Kc = 0.2  # Dissociation constant for calcium effect on release probability
m = 4  # Hill coefficient for calcium effect on release probability
tau_f = 0.1
tau_d = 0.6  # Increased tau_d to reduce depression
tau_DA = 0.5
tau_rec = 0.3  # Reduced tau_rec to speed up vesicle replenishment
tau_ipsc = 0.02  # Time constant for IPSC decay
dt = 0.001
time = np.arange(0, 2, dt)
n_infty = 1.0  # Steady-state vesicle availability
spike_times = [0.5, 0.6]  # Example spike times in seconds (short interval for paired-pulse)

# New parameters for calcium dynamics
calcium_influx = 0.15  # Increased calcium influx
calcium_channel_open_prob = 0.1  # Increased calcium channel open probability

# Observed PPR values for each condition
PPR_data_no_agonist = {
    'control': 0.9,
    '6-OHDA': 1
}

PPR_data_with_agonist = {
    'Quinpirole': 10,
    'Quinpirole 6-OHDA': 0.1
}

# Dopamine concentrations for each condition
DA_concentrations = {
    'control': 30e-9,
    '6-OHDA': 0.0,
    'Quinpirole':30e-9,
    'Quinpirole 6-OHDA': 0.0
}

# Receptor effects for each condition
D2_receptor_effects = {
    'control': 1.0,
    '6-OHDA': 1.0,
    'Quinpirole': 100000,  # Increased effect
    'Quinpirole 6-OHDA': 100000  # Increased effect
}

def synaptic_response(params, condition):
    Ca0, Kc, tau_f, tau_d, tau_rec,P0, calcium_influx, calcium_channel_open_prob = params
    DA_concentration = DA_concentrations[condition]
    D2_effect = D2_receptor_effects[condition]
    
    u = np.zeros_like(time)
    R = np.ones_like(time)
    Ca = np.zeros_like(time)
    P = np.zeros_like(time)
    A = np.zeros_like(time)
    n = np.ones_like(time)  # Vesicle availability
    IPSC = np.zeros_like(time)

    u[0] = P0
    R[0] = 1
    n[0] = n_infty

    for t in range(1, len(time)):
        Ca[t] = Ca0 + (DA_concentration**n_Hill / (Kd_DA**n_Hill + DA_concentration**n_Hill)) * D2_effect + calcium_influx
        P[t] = P0 * (Ca[t]**m / (Kc**m + Ca[t]**m)) * calcium_channel_open_prob 
        dn = (n_infty - n[t-1]) / tau_rec * dt
        if time[t] in spike_times:
            dn -= u[t-1] * R[t-1] * n[t-1]
        n[t] = n[t-1] + dn
        
        du = (P[t] - u[t-1]) / tau_f * dt
        if time[t] in spike_times:
            du += P[t] * (1 - u[t-1])
        u[t] = u[t-1] + du
        
        dR = (1 - R[t-1]) / tau_d * dt
        if time[t] in spike_times:
            dR -= u[t-1] * R[t-1] * n[t-1]
        R[t] = R[t-1] + dR

        # Clamp values to prevent overflow
        u[t] = np.clip(u[t], 0, 1)
        R[t] = np.clip(R[t], 0, 1)
        n[t] = np.clip(n[t], 0, n_infty)
        
        if time[t] in spike_times:
            A[t] = u[t] * R[t] * n[t]
            IPSC[t:] += A[t] * np.exp(-(time[t:] - time[t]) / tau_ipsc)

    A1 = A[time == spike_times[0]]
    A2 = A[time == spike_times[1]]
    
    PPR_model = A2 / A1 if A1 != 0 else np.inf
    
    return PPR_model if isinstance(PPR_model, np.ndarray) else np.array([PPR_model]), IPSC

def objective_function_no_agonist(params):
    total_error = 0
    for condition, observed_PPR in PPR_data_no_agonist.items():
        PPR_model, _ = synaptic_response(params, condition)
        total_error += (PPR_model[0] - observed_PPR) ** 2
    return total_error

def objective_function_with_agonist(modified_params, fixed_params):
    total_error = 0
    full_params = fixed_params.copy()
    full_params = modified_params  # Modify specific parameters

    for condition, observed_PPR in PPR_data_with_agonist.items():
        PPR_model, _ = synaptic_response(full_params, condition)
        total_error += (PPR_model[0] - observed_PPR) ** 2
    return total_error

# Initial guess for parameters
initial_params_no_agonist = [Ca0, Kc, tau_f, tau_d, tau_rec, calcium_influx, P0, calcium_channel_open_prob]

# Define bounds for each parameter to ensure they remain within physiological ranges
bounds_no_agonist = [
    (1e-9, 1e-5),  # Ca0: cannot be negative
    (0, None),  # Kc: cannot be negative
    (0, None),  # tau_rec: cannot be negative
    (0.01, 1),  # P0: probability between 0 and 1
    (0.01, None),  # calcium_influx: cannot be negative
    (0, 0.8),  # calcium_channel_open_prob: probability between 0 and 1
    (0, None),  # tau_f: cannot be negative
    (0, None),  # tau_d: cannot be negative
]

bounds_with_agonist = [
    (0, 1),  # P0: probability between 0 and 1
    (0, None),  # calcium_influx: cannot be negative
    (0, 1) , # calcium_channel_open_prob: probability between 0 and 1
    (0, None),  # tau_f: cannot be negative
    (0, None)  # tau_d: cannot be negative
]

# Optimize parameters to fit the observed PPR for conditions without agonist
result_no_agonist = minimize(objective_function_no_agonist, initial_params_no_agonist, method='L-BFGS-B', bounds=bounds_no_agonist)
fitted_params_no_agonist = result_no_agonist.x
initial_modified_params = fitted_params_no_agonist

# Print the fitted parameters for conditions without agonist
print("Fitted Parameters without agonist:", fitted_params_no_agonist)

# Optimize parameters to fit the observed PPR for conditions with agonist by modifying specific parameters
result_with_agonist = minimize(objective_function_with_agonist, initial_params_no_agonist, args=(fitted_params_no_agonist,), method='L-BFGS-B', bounds=bounds_no_agonist)
fitted_params_with_agonist = result_with_agonist.x

# Combine the optimized parameters
fitted_params_combined = fitted_params_no_agonist.copy()
fitted_params_combined = fitted_params_with_agonist

# Print the fitted parameters for conditions with agonist
print("Fitted Parameters with agonist:", fitted_params_with_agonist)

# Plot the model responses with fitted parameters for all conditions
conditions = ['control', '6-OHDA', 'Quinpirole', 'Quinpirole 6-OHDA']

conditions_no_agonist = ['control', '6-OHDA']
plt.figure(figsize=(12, 6))

for i, condition in enumerate(conditions_no_agonist):
    _, IPSC = synaptic_response(fitted_params_no_agonist, condition)
    plt.subplot(2, 1, i+1)
    plt.plot(time, IPSC, label=f'IPSC - {condition}')
    plt.legend()

plt.tight_layout()
plt.show()



# Plot the model responses for conditions with agonist
conditions_with_agonist = ['Quinpirole', 'Quinpirole 6-OHDA']
plt.figure(figsize=(12, 6))

for i, condition in enumerate(conditions_with_agonist):
    _, IPSC = synaptic_response(fitted_params_combined, condition)
    plt.subplot(2, 1, i+1)
    plt.plot(time, IPSC, label=f'IPSC - {condition}')
    plt.legend()

plt.tight_layout()
plt.show()