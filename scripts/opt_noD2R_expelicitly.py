import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution, basinhopping
from scipy.stats import mannwhitneyu 
import utils 
import warnings
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in scalar divide")

data_file="results/features_all_data.pkl"
all_data=pd.read_pickle(data_file)
# Define the model parameters and time constants
Kd_DA = 10e-9  # Dissociation constant for D2 receptor (10 nM converted to M)
n_Hill = 4  # Hill coefficient for DA effect on calcium
P0 = 0.1
Kc = 1e-6   # Dissociation constant for calcium effect on release probability
m = 4  # Hill coefficient for calcium effect on release probability
tau_f = 0.2
tau_d = 0.9  # Increased tau_d to reduce depression
tau_rec = 0.5  # Reduced tau_rec to speed up vesicle replenishment
tau_ipsc = 0.02  # Time constant for IPSC decay
dt = 0.001
time = np.arange(0, 1, dt)
n_infty = 1.0  # Steady-state vesicle availability
spike_times = [0.3, 0.4]  # Example spike times in seconds (short interval for paired-pulse)
Ca0=50e-9
# New parameters for calcium dynamics
calcium_influx =5e-6    # Increased calcium influx
calcium_channel_open_prob = 0.2  # Increased calcium channel open probability
tau_Ca = 0.1  # Calcium decay time constant (in seconds)

# Observed PPR values for each condition
PPR_data = {
    'control': 1,
    '6-OHDA': 0.8,
    'Quinpirole': 2,
    'Quinpirole 6-OHDA': 1.
}

# Dopamine concentrations for each condition
DA_concentrations = {
    'control': 50e-9,
    '6-OHDA': 1e-11,
    'Quinpirole': 50e-9,
    'Quinpirole 6-OHDA': 0.0
}

quinpirole_concentrations = {
    'control': 0.0,
    '6-OHDA': 0.0,
    'Quinpirole': 10e-6,           # 10 μM quinpirole
    'Quinpirole 6-OHDA': 10e-6
}


def synaptic_response(params, condition):
    P0, Kc, tau_f, tau_d, tau_rec, tau_Ca, calcium_influx, calcium_channel_open_prob = params
    DA_concentration = DA_concentrations[condition]
    agonist_concentration = quinpirole_concentrations[condition]
    D2R_activation = utils.D2R_activation(DA_concentration, agonist_concentration)
    # Use the defined D2R_activation function

    u = np.zeros_like(time)
    R = np.ones_like(time)
    Ca = np.zeros_like(time)
    P = np.zeros_like(time)
    A = np.zeros_like(time)
    n = np.ones_like(time)  # Vesicle availability
    IPSC = np.zeros_like(time)
    Ca[0] = Ca0
    D2R_activation_trace = np.zeros_like(time)

    u[0] = P0
    R[0] = 1
    n[0] = n_infty

    for t in range(1, len(time)):
        dCa_dt = - (Ca[t-1] - Ca0) / tau_Ca
        if utils.is_spike_time(time[t], spike_times, dt):
            inhibition_factor = 0 
            inhibition_factor = np.clip(inhibition_factor, 0, 1)
            Ca_influx_effective = calcium_influx 
            dCa_dt += Ca_influx_effective / tau_Ca  # Adjust units appropriately

        Ca[t] = Ca[t-1] + dCa_dt * dt
        P[t] = P0 * (Ca[t]**m / (Kc**m + Ca[t]**m)) * calcium_channel_open_prob 
        P[t] = np.clip(P[t], 0, 1)  # Clamp P(t) to avoid extreme values
        dn = (n_infty - n[t-1]) / tau_rec * dt

        if time[t] in spike_times:
            dn -= u[t-1] * R[t-1] * n[t-1]
        n[t] = n[t-1] + dn
        
        du = (- u[t-1]) / tau_f * dt
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
            IPSC[t:] +=  A[t] * np.exp(-(time[t:] - time[t]) / tau_ipsc)
    
    D2R_activation_trace[t] = D2R_activation

    idx1 = np.where(np.isclose(time, spike_times[0], atol=dt/2))[0]
    idx2 = np.where(np.isclose(time, spike_times[1], atol=dt/2))[0]

    if idx1.size > 0 and idx2.size > 0:
        A1 = A[idx1[0]]
        A2 = A[idx2[0]]
        if A1 == 0:
            PPR_model = np.inf
        else:
            PPR_model = A2 / A1
    else:
        PPR_model = np.nan  # Handle cases where spikes are not found

    # Print intermediate values
    #print(f"Condition: {condition}, PPR_model: {PPR_model}, A1: {A1}, A2: {A2}")

    
    return np.array([PPR_model]), IPSC, Ca, D2R_activation_trace,  np.trapz(P, time),  np.trapz(Ca, time)

def objective_function(params, condition='control'):
    total_error = 0
    observed_PPR=PPR_data[condition]
    PPR_model, IPSC_model,_,_,_,_ = synaptic_response(params, condition)
    total_error = (PPR_model[0] - observed_PPR) ** 2
    IPSC_error = (np.max(IPSC_model) - observed_IPSC_quinpirole) ** 2
    # Calculate the peak IPSC and add penalty if it exceeds 2500 pA

    return total_error 

def objective_function_q(params, condition='Quinpirole'):
    total_error = 0
    observed_PPR=PPR_data[condition]
    PPR_model, IPSC_model,_,_,_,_ = synaptic_response(params, condition)
    total_error = (PPR_model[0] - observed_PPR) ** 2
    IPSC_error = (np.max(IPSC_model) - observed_IPSC_quinpirole) ** 2
    # Calculate the peak IPSC and add penalty if it exceeds 2500 pA

    return total_error 

def objective_function_quinpirole(params):
    # Unpack parameters to optimize
    P0, calcium_influx, tau_f, tau_d, tau_Ca, calcium_channel_open_prob = params

    # Fixed parameters from control condition
    Kc = control_params_dict['Kc']
    tau_rec = control_params_dict['tau_rec']
    m = 4  # Assuming m is defined globally

    # Model parameters for quinpirole condition
    model_params_quinpirole = [P0, Kc, tau_f, tau_d, tau_rec, tau_Ca, calcium_influx, calcium_channel_open_prob]

    # Run the model for quinpirole condition
    PPR_model_quinpirole, IPSC_quinpirole, Ca_quinpirole,_,_,_= synaptic_response(model_params_quinpirole, 'Quinpirole')

    # Run the model for control condition using fixed parameters
    model_params_control = [
        control_params_dict['P0'],
        control_params_dict['Kc'],
        control_params_dict['tau_f'],
        control_params_dict['tau_d'],
        control_params_dict['tau_rec'],
        control_params_dict['tau_Ca'],
        control_params_dict['calcium_influx'],
        control_params_dict['calcium_channel_open_prob']
    ]
    PPR_model_control, IPSC_control, Ca_control,_,_,_ = synaptic_response(model_params_control, 'control')

    # Extract the peak IPSC amplitudes
    peak_IPSC_quinpirole = np.max(IPSC_quinpirole)
    peak_IPSC_control = np.max(IPSC_control)

    # Get observed PPR values
    observed_PPR_quinpirole = PPR_data['Quinpirole']

    # Calculate error between model and experimental data
    PPR_error = (PPR_model_quinpirole[0] - observed_PPR_quinpirole) ** 2

    # Initialize amplitude penalty
    amplitude_penalty = 0
    amplitude_penalty_= 0
    # First Condition: Penalty if quinpirole peak is greater than control peak
    if peak_IPSC_quinpirole > peak_IPSC_control:
        amplitude_penalty += (peak_IPSC_quinpirole - peak_IPSC_control) ** 2 * 1e3  # Adjust penalty weight as needed

    # Second Condition: Penalty if control peak is more than 3 times larger than quinpirole peak
    if peak_IPSC_quinpirole != 0:
        ratio = peak_IPSC_control / peak_IPSC_quinpirole
    else:
        ratio = np.inf  # Assign a large value to apply penalty

    if ratio > 3:
         amplitude_penalty_ += ratio ** 5 * 1e3  # Adjust penalty weight as needed

    total_error = PPR_error + amplitude_penalty_

    return total_error
# Initial guess for parameters
initial_params = [P0, Kc, tau_f, tau_d, tau_rec,tau_Ca, calcium_influx, calcium_channel_open_prob]
initial_params_q = [P0, Kc, tau_f, tau_d, tau_rec,tau_Ca, calcium_influx, calcium_channel_open_prob]

# Define bounds for each parameter to ensure they remain within physiological ranges
bounds = [
    (0.1   , 1),  # P0: probability between 0 and 1
    (1e-9, 1e-5), # Kc: between 1 nM and 10 µM (adjust as needed)
    (0.00, 1),  # tau_f: cannot be negative
    (0.00, 2),  # tau_d: cannot be negative
    (0.00, 2),  # tau_rec: cannot be negative
    (2e-2, 5e-1),  # tau_ca: calcium decay 
    (1e-7, 100e-5),  # calcium_influx: cannot be negative
    (0, 0.6)  # calcium_channel_open_prob: probability between 0 and 1
]


# Initialize a dictionary to store IPSC traces
IPSC_data = {}
Conditions=["10uM QP","6-OHDA","10uM CdCl2"]
#Regions=["DL","VM","DM","VL"]
Regions=["VM"]

Phase=["Baseline","During"]

fitted_params_control_list = []
fitted_params_quinpirole_list = []
cell_indices = []

# Ensure that ppr_values and peak_amp_values are aligned
# For simplicity, assume they are in the same order and correspond to each cell$
fitted_params_control_list = []
fitted_params_quinpirole_list = []

for region in Regions:
    ppr_values = all_data["PPR"].loc["10uM QP"].loc[region]["Baseline"].dropna()
    peak_amp_values = all_data["Peak amplitude"].loc["10uM QP"].loc[region]["Baseline"].dropna()
    peak_amp_quinp = all_data["Peak amplitude"].loc["10uM QP"].loc[region]["During"].dropna()
    ppr_values_quinpirole = all_data["PPR"].loc["10uM QP"].loc[region]["During"].dropna()
    
    control_iterations = []
    control_params_history = []
    print(region)
    print(ppr_values)
    for idx, (ppr, ppr_quinp,peak_amp,peak_amp_q) in enumerate(zip(ppr_values, ppr_values_quinpirole,peak_amp_values,peak_amp_quinp)):
        observed_PPR = ppr
        observed_PPR_quinpirole = ppr_quinp  # Assuming PPR under quinpirole is the same as observed
        observed_IPSC = peak_amp  # Use peak amplitude as observed IPSC
        observed_IPSC_quinpirole = peak_amp_q  # Use peak amplitude as observed IPSC


        # Update PPR_data for control and quinpirole
        PPR_data['control'] = observed_PPR  # Assuming you have control PPR data
        PPR_data['Quinpirole'] = observed_PPR_quinpirole

        if observed_PPR > 1 :
            initial_params[0]=0.2 + np.random.uniform(0, 0.2)#P_0=
            initial_params[2]=1+ np.random.uniform(-0.2, 0.2)    # tau_f
            initial_params[3]=2+ np.random.uniform(-0.5, 0.5)    #tau_D
            #bounds[0]=(0.1, 0.5)
        else: 
            initial_params[0]=0.8+ np.random.uniform(-0.2, 0.1)   
            initial_params[2]=0.05+ np.random.uniform(-0.05, 0.05)   
            initial_params[3]=0.1+ np.random.uniform(-0.05, 0.1)
            #bounds[0]=(0.5, 1)

        # Perform optimization
        result_control = minimize(objective_function, initial_params, method='L-BFGS-B', bounds=bounds)
        
        control_params = result_control.x

        PPR_model_cont, IPSC_cont, Ca_cont, D2R_activation_cont, tot_P_c, tot_Ca_c = synaptic_response(control_params, 'control')
        #if np.max(Ca_cont) > 1e-6:
        #    continue  # Skip the rest of this loop iteration
        control_params_dict = {
            'Region': region,
            'P0': control_params[0],
            'Kc': control_params[1],
            'tau_f': control_params[2],
            'tau_d': control_params[3],
            'tau_rec': control_params[4],
            'tau_Ca': control_params[5],
            'calcium_influx': control_params[6],
            'calcium_channel_open_prob': control_params[7],
            'total_P': tot_P_c,  # Total release probability for quinpirole
            'total_Ca': tot_Ca_c 
        }
        fitted_params_control_list.append(control_params_dict)

        if observed_PPR_quinpirole > 1.05 :

            initial_params_q[0]=0.2 + np.random.uniform(0, 0.2)#P_0=
            initial_params_q[2]=1+ np.random.uniform(-0.2, 0.2)    # tau_f
            initial_params_q[3]=2+ np.random.uniform(-0.5, 0.5)    #tau_D
            bounds_quinpirole[0]=(0.1, 0.5)

        else: 

            initial_params_q[0]=0.8+ np.random.uniform(-0.2, 0.1)   
            initial_params_q[2]=0.05+ np.random.uniform(-0.05, 0.05)   
            initial_params_q[3]=0.1+ np.random.uniform(-0.05, 0.1)
            bounds_quinpirole[0]=(0.5, 1)

        # Optimize for quinpirole condition
        result_quinpirole = minimize(
            objective_function_q,
            initial_params_q,
            method='L-BFGS-B',
            bounds=bounds)
        
        differential_evolution(objective_function, bounds_quinpirole, strategy='best1bin', maxiter=100)  
        fitted_params_quinpirole = result_quinpirole.x

        # Run the synaptic response to get tot_P and tot_Ca for quinpirole
        PPR_model_quinpirole, IPSC_quinpirole, Ca_quinpirole, D2R_activation_quinpirole, tot_P_q, tot_Ca_q = synaptic_response(fitted_params_quinpirole, 'Quinpirole')

        fitted_params_quinpirole_dict = {
            'Region': region,
            'P0': fitted_params_quinpirole[0],
            'Kc': fitted_params_quinpirole[1],
            'tau_f': fitted_params_quinpirole[2],
            'tau_d': fitted_params_quinpirole[3],
            'tau_rec': fitted_params_quinpirole[4],
            'tau_Ca': fitted_params_quinpirole[5],
            'calcium_influx': fitted_params_quinpirole[6],
            'calcium_channel_open_prob': fitted_params_quinpirole[7],
            'total_P': tot_P_q,  # Total release probability for quinpirole
            'total_Ca': tot_Ca_q 
        }
        fitted_params_quinpirole_list.append(fitted_params_quinpirole_dict)

        # Store the cell index
        cell_indices.append(idx)

print(cell_indices)

for idx in cell_indices:
    # Get the control and quinpirole parameters
    control_params_dict = fitted_params_control_list[idx]
    fitted_params_quinpirole_dict = fitted_params_quinpirole_list[idx]

    # Reconstruct full parameter lists
    model_params_control = [
        control_params_dict['P0'],
        control_params_dict['Kc'],
        control_params_dict['tau_f'],
        control_params_dict['tau_d'],
        control_params_dict['tau_rec'],
        control_params_dict['tau_Ca'],
        control_params_dict['calcium_influx'],
        control_params_dict['calcium_channel_open_prob']
    ]

    model_params_quinpirole = [
        fitted_params_quinpirole_dict['P0'],
        fitted_params_quinpirole_dict['Kc'],  # Kc remains the same
        fitted_params_quinpirole_dict['tau_f'],
        fitted_params_quinpirole_dict['tau_d'],
        fitted_params_quinpirole_dict['tau_rec'],  # tau_rec remains the same
        fitted_params_quinpirole_dict['tau_Ca'],  # tau_Ca remains the same
        fitted_params_quinpirole_dict['calcium_influx'],
        fitted_params_quinpirole_dict['calcium_channel_open_prob']
    ]

    # Run the model for control and quinpirole
    PPR_model_control, IPSC_control , Ca_control, D2R_activation_control,_,_ = synaptic_response(model_params_control, 'control')
    PPR_model_quinpirole, IPSC_quinpirole,  Ca_quinpirole, D2R_activation_quinpirole,_,_ = synaptic_response(model_params_quinpirole, 'Quinpirole')
    # Get the observed PPR values
    observed_PPR_control = ppr_values.iloc[idx]
    observed_PPR_quinpirole = ppr_values_quinpirole.iloc[idx]

    # Plot the synaptic responses
    # Plot the synaptic responses and traces
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # IPSC plots
    axes[0, 0].plot(time, IPSC_control, label='Model IPSC')
    axes[0, 0].set_title(f'Cell {idx + 1} - Control Condition')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('IPSC')
    axes[0, 0].legend()
    axes[0, 0].text(0.05, 0.95, f"Observed PPR: {observed_PPR_control:.2f}\nSimulated PPR: {PPR_model_control[0]:.2f} \nPr: {control_params_dict['P0']:.2f}",
                    transform=axes[0, 0].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    axes[0, 1].plot(time, IPSC_quinpirole, label='Model IPSC', color='orange')
    axes[0, 1].set_title(f'Cell {idx + 1} - Quinpirole Condition')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('IPSC')
    axes[0, 1].legend()
    axes[0, 1].text(0.05, 0.95, f"Observed PPR: {observed_PPR_quinpirole:.2f}\nSimulated PPR: {PPR_model_quinpirole[0]:.2f}\nPr: {fitted_params_quinpirole_dict['P0']:.2f}",
                    transform=axes[0, 1].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # Calcium traces
    axes[1, 0].plot(time, Ca_control, label='Calcium Trace')
    axes[1, 0].set_title('Calcium Concentration - Control')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Calcium Concentration (M)')
    axes[1, 0].legend()

    axes[1, 1].plot(time, Ca_quinpirole, label='Calcium Trace', color='orange')
    axes[1, 1].set_title('Calcium Concentration - Quinpirole')
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Calcium Concentration (M)')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.show()



df_control = pd.DataFrame(fitted_params_control_list)
df_control['Condition'] = 'Control'
df_quinpirole = pd.DataFrame(fitted_params_quinpirole_list)
df_quinpirole['Condition'] = 'Quinpirole'

df_all = pd.concat([df_control, df_quinpirole], ignore_index=True)

# List of parameters to analyze
parameters = ['calcium_channel_open_prob', 'calcium_influx', 'tau_Ca', 'total_Ca', 'P0']

# Loop over each region and plot
for region in Regions:
    fig, axes = plt.subplots(len(parameters), 1, figsize=(10, len(parameters)*4))
    region_data = df_all

    for i, param in enumerate(parameters):
        sns.violinplot(data=region_data, x='Condition', y=param, ax=axes[i], inner=None, palette="Set2")
        sns.scatterplot(data=region_data, x='Condition', y=param, ax=axes[i], color="black", s=50)
        axes[i].set_title(f'{param} - {region}')

        # Statistical test (Mann-Whitney U test)
        control_data = region_data[region_data['Condition'] == 'Control'][param]
        quinpirole_data = region_data[region_data['Condition'] == 'Quinpirole'][param]
        stat, p_value = mannwhitneyu(control_data, quinpirole_data)

        # Mark statistical significance
        if p_value < 0.001:
            significance = '***'
        elif p_value < 0.01:
            significance = '**'
        elif p_value < 0.05:
            significance = '*'
        else:
            significance = 'ns'  # not significant

        axes[i].text(0.5, 1.05, f'p = {p_value:.4f} {significance}', ha='center', va='bottom', 
                     transform=axes[i].transAxes, fontsize=12)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"Stat_{region}_old.pdf")