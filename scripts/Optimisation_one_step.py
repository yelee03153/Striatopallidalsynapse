import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution, basinhopping
from scipy.stats import mannwhitneyu, wilcoxon
import utils 
import warnings
import pandas as pd
import seaborn as sns

def callback_control(xk, convergence):
    iteration = len(control_iterations)
    control_iterations.append(iteration)
    control_params_history.append(xk.copy())

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
calcium_channel_open_prob = 0.3  # Increased calcium channel open probability
tau_Ca = 0.1  # Calcium decay time constant (in seconds)

PPR_data = {
    'control': 1,
    '6-OHDA': 0.8,
    'Quinpirole': 2,
    'Quinpirole 6-OHDA': 1.
}# Dopamine concentrations for each condition
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


def synaptic_response(params, condition,scaling_params = None):
    P0, Kc, tau_f, tau_d, tau_rec, tau_Ca, calcium_influx, calcium_channel_open_prob = params
    #DA_concentration = DA_concentrations[condition]
    #agonist_concentration = quinpirole_concentrations[condition]
    #D2R_activation = utils.D2R_activation(DA_concentration, agonist_concentration)
    #inhibition_factor = np.clip(D2R_activation, 0, 1)
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

    total_release_probability = 0.0


    if scaling_params is not None and condition == 'Quinpirole':
        scale_open_prob, scale_d, scale_f= scaling_params
        calcium_channel_open_prob *=  scale_open_prob
        #calcium_influx *=  scale_influx
        #tau_d *= scale_d
        #tau_f *= scale_f

    for t in range(1, len(time)):
        dCa_dt = - (Ca[t-1] - Ca0) / tau_Ca
        if utils.is_spike_time(time[t], spike_times, dt):
            Ca_influx_effective = calcium_influx* calcium_channel_open_prob

            dCa_dt += Ca_influx_effective / tau_Ca  # Adjust units appropriately

        Ca[t] = Ca[t-1] + dCa_dt * dt
        P[t] = P0 * (Ca[t]**m / (Kc**m + Ca[t]**m)) 
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
    
    #D2R_activation_trace[t] = D2R_activation

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
        PPR_model = np.nan  # HanVLe cases where spikes are not found

    # Print intermediate values
    #print(f"Condition: {condition}, PPR_model: {PPR_model}, A1: {A1}, A2: {A2}")

    
    return np.array([PPR_model]), IPSC, P, Ca,  np.trapz(P, time),  np.trapz(Ca, time)

def objective_function(params, condition='control'):
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

    if ratio > 4:
         amplitude_penalty_ += ratio ** 5 * 1e3  # Adjust penalty weight as needed

    total_error = PPR_error + amplitude_penalty_

    return total_error

def objective_function_quinpirole_p(scaling_params):
    observed_PPR_quinpirole=PPR_data['Quinpirole']
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
    PPR_model_quinpirole, IPSC_quinpirole, Ca_quinpirole,_,_,_= synaptic_response(model_params_control, 'Quinpirole',scaling_params =scaling_params)
    error = (PPR_model_quinpirole[0] - observed_PPR_quinpirole) ** 2

    _, IPSC_control, _,_,_,_ = synaptic_response(model_params_control, 'control')

    # Extract the peak IPSC amplitudes
    peak_IPSC_quinpirole = np.max(IPSC_quinpirole)
    peak_IPSC_control = np.max(IPSC_control)
    ratio = peak_IPSC_control / peak_IPSC_quinpirole
    amplitude_penalty_=0
    if ratio > 4:
         amplitude_penalty_ += ratio ** 5 * 1e3 

    return error +amplitude_penalty_


# Initial guess for parameters
initial_params = [P0, Kc, tau_f, tau_d, tau_rec,tau_Ca, calcium_influx, calcium_channel_open_prob]

# Define bounds for each parameter to ensure they remain within physiological ranges
bounds = [
    (0.1   , 0.9),  # P0: probability between 0 and 1
    (1e-9, 1e-5), # Kc: between 1 nM and 10 µM (adjust as needed)
    (0.00, 1),  # tau_f: cannot be negative
    (0.00, 2),  # tau_d: cannot be negative
    (0.00, 2),  # tau_rec: cannot be negativels -
    (2e-2, 5e-1),  # tau_ca: calcium decay 
    (1e-7, 100e-5),  # calcium_influx: cannot be negative
    (0.1   , 0.7)  # calcium_channel_open_prob: probability between 0 and 1
]


# Initialize a dictionary to store IPSC traces
IPSC_data = {}
Conditions=["10uM QP","6-OHDA","10uM CdCl2"]
#Regions=["VL","VL","VL","VL"]
Regions=["VL"]

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
        observed_PPR = np.mean(ppr_values)
        observed_PPR_quinpirole = np.mean(ppr_values_quinpirole)  # Assuming PPR under quinpirole is the same as observed
        observed_IPSC = peak_amp  # Use peak amplitude as observed IPSC
        observed_IPSC_quinpirole = peak_amp_q  # Use peak amplitude as observed IPSC


        # Update PPR_data for control and quinpirole
        PPR_data['control'] = observed_PPR  # Assuming you have control PPR data
        PPR_data['Quinpirole'] = observed_PPR_quinpirole

        # Optimize for control condition


        #initial_params[7]= 0.2 + np.random.uniform(-0.1, 0.1) 

        control_params = None
        best_control_score = float('inf')

        # Multi-start optimization for control condition
        #for start in range(10):
             #Create random initial parameters based on some reasonable base
        """
        initial_params = [
            0.5 + np.random.uniform(-0.1, 0.1),  # P0
            1e-6,  # Kc
            0.5 + np.random.uniform(-0.4, 0.4),  # tau_f
            1.0 + np.random.uniform(-0.5, 0.5),  # tau_d
            0.5,  # tau_rec
            0.1 + np.random.uniform(-0.01, 0.01),  # tau_Ca
            5e-6 + np.random.uniform(-1e-6, 1e-6),  # calcium_influx
            0.2 + np.random.uniform(0, 0.2)  # calcium_channel_open_prob
        ]
        """
        if observed_PPR > 1.05 :
            initial_params[0]=0.2 + np.random.uniform(0, 0.2)#P_0=
            initial_params[2]=1+ np.random.uniform(-0.2, 0.2)    # tau_f
            initial_params[3]=2+ np.random.uniform(-0.5, 0.5)    #tau_D
            bounds[0]=(0.1, 0.5)
        else: 
            initial_params[0]=0.8+ np.random.uniform(-0.2, 0.1)   
            initial_params[2]=0.05+ np.random.uniform(-0.05, 0.05)   
            initial_params[3]=0.1+ np.random.uniform(-0.05, 0.1)
            bounds[0]=(0.55, 0.9)

            # Perform optimization
        #result_control = minimize(objective_function, initial_params, method='L-BFGS-B', bounds=bounds)
            
            # Check if this is the best score
            #if result_control.fun < best_control_score:
            #    best_control_score = result_control.fun
            #    control_params = result_control.x

        print(bounds[0])

        result_control = differential_evolution(objective_function, bounds, strategy='best1bin', maxiter=10)  
        #result_control = minimize(objective_function, initial_params, method='L-BFGS-B', bounds=bounds)
        #result_control =basinhopping(objective_function, initial_params, niter=50)
        control_params = result_control.x
                # Perform global optimization

        # Use the result from the global optimizer as initial parameters for the local optimizer


        PPR_model_cont, IPSC_cont, Ca_cont, D2R_activation_cont, tot_P_c, tot_Ca_c = synaptic_response(control_params, 'control')
        print(observed_PPR)

        print(PPR_model_cont)

        #if np.max(Ca_cont) < 1e-6:
        #    continue  # Skip the rest of this loop iteration

        # Set initial parameters and bounds for quinpirole optimization
        """
        initial_params_quinpirole = [
            P0,
            control_params_dict['calcium_influx'],
            tau_f,
            tau_d,
            control_params_dict['tau_Ca'],
            0.2
        ]
        """
        bounds_quinpirole = [
            (0.1, 0.9),
            (1e-7, 100e-5),
            (0.001, 1.5),
            (0.001, 2.5),
            (2e-2, 5e-1),
            (0, control_params[7]+0.05)  # calcium_channel_open_prob: probability between 0 and 1
        ]

        initial_scaling_params = [0.2,  1, 1]  # Start at 1 (neutral scaling)
        bounds_scaling = [(0.01, 1.05), (0.8, 1.2), (0.8, 1.2)]  # Limits for scaling factors
        if observed_PPR < 1.08 and observed_PPR_quinpirole > 1.05:
            print(control_params[2])
            print(control_params[3])

        # Synapse changes from facilitating to depressing; P0 should increase
            control_params[0]=0.2
            control_params[0]=np.clip(control_params[0], 0.1, 0.5)
            bounds_scaling[1]=(0.1,55)
            bounds_scaling[2]=(0.1, 55)
            control_params[2] =0.2

            control_params[3] = 1

            print(control_params[0])


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


        """
        if observed_PPR_quinpirole > 1.05 :
            initial_params_quinpirole[2]=1 +np.random.uniform(-0.2, 0.2) # tau_f
            initial_params_quinpirole[3]=2+ np.random.uniform(-0.5, 0.5)   #tau_D
            initial_params_quinpirole[0]=0.2+ np.random.uniform(0, 0.2)
            bounds_quinpirole[0]=(0.1, 0.5)

        else: 
            initial_params_quinpirole[2]=0.05 + np.random.uniform(-0.05, 0.05)   
            initial_params_quinpirole[3]=0.1+ np.random.uniform(-0.05, 0.1)   
            initial_params_quinpirole[0]=0.8+ np.random.uniform(-0.2, 0.1)
            bounds_quinpirole[0]=(0.55, 0.9)
        """ 
        #check if syn cahnges modality under quinperole 

        # Optimize for quinpirole condition
        #result_quinpirole = minimize(
        #    objective_function_quinpirole,
        #    initial_params_quinpirole,
        #    method='L-BFGS-B',
        #    bounds=bounds_quinpirole)
        

        result_quinpirole = differential_evolution(objective_function_quinpirole_p, bounds_scaling, strategy='best1bin', maxiter=300) 
        #result_quinpirole = minimize(objective_function_quinpirole_p, initial_scaling_params,  method='L-BFGS-B',  bounds=bounds_scaling) 

        #fitted_params_quinpirole = result_quinpirole.x
        scaling_params_quinpirole = result_quinpirole.x

        """
        model_params_quinpirole = [
            fitted_params_quinpirole[0],  # P0 (optimized)
            control_params_dict['Kc'],  # Kc (fixed from control)
            fitted_params_quinpirole[2],  # tau_f (optimized)
            fitted_params_quinpirole[3],  # tau_d (optimized)
            control_params_dict['tau_rec'],  # tau_rec (fixed from control)
            fitted_params_quinpirole[4],  # tau_Ca (optimized)
            fitted_params_quinpirole[1],  # calcium_influx (optimized)
            fitted_params_quinpirole[5]  # calcium_channel_open_prob (optimized)
        ]

        """

        # Run the synaptic response to get tot_P and tot_Ca for quinpirole
        PPR_model_quinpirole, IPSC_quinpirole, Ca_quinpirole, _, tot_P_q, tot_Ca_q = synaptic_response(control_params, 'Quinpirole', scaling_params=scaling_params_quinpirole)
        print(observed_PPR_quinpirole)

        print(PPR_model_quinpirole)
        
        if (PPR_model_quinpirole[0] > 1 and observed_PPR_quinpirole < 1) or (PPR_model_quinpirole[0] < 1 and observed_PPR_quinpirole > 1) or abs(PPR_model_quinpirole-observed_PPR_quinpirole)>0.4:
            print(f"Excluding cell {idx} due to mismatch in PPR behavior.")
            continue  # Skip adding this cell to both control and quinpirole lists

        fitted_params_quinpirole_dict = {
            'Region': region,
            'P0': control_params[0],
            "Kc": control_params[1],
            'calcium_influx': control_params[6],
            'tau_f': control_params[2],
            'tau_d': control_params[3],
            'tau_rec': control_params[4],
            'tau_Ca': control_params[5],
            'calcium_channel_open_prob': control_params[7]*scaling_params_quinpirole[0],
            'total_P': tot_P_q,  # Total release probability for quinpirole
            'total_Ca': tot_Ca_q  # Total calcium for quinpirole
        }
        fitted_params_quinpirole_list.append(fitted_params_quinpirole_dict)
        fitted_params_control_list.append(control_params_dict)

        # Store the cell index
        cell_indices.append(idx)



for i,idx in enumerate(cell_indices):
    # Get the control and quinpirole parameters
    control_params_dict = fitted_params_control_list[i]
    fitted_params_quinpirole_dict = fitted_params_quinpirole_list[i]

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
    """
    model_params_quinpirole = [
        fitted_params_quinpirole_dict['P0'],
        control_params_dict['Kc'],  # Kc remains the same
        fitted_params_quinpirole_dict['tau_f'],
        fitted_params_quinpirole_dict['tau_d'],
        control_params_dict['tau_rec'],  # tau_rec remains the same
        fitted_params_quinpirole_dict['tau_Ca'],  # tau_Ca remains the same
        fitted_params_quinpirole_dict['calcium_influx'],
        fitted_params_quinpirole_dict['calcium_channel_open_prob']
    ]
    """
    model_params_quinpirole = [
        fitted_params_quinpirole_dict['P0'],
        fitted_params_quinpirole_dict['Kc'],
        fitted_params_quinpirole_dict['tau_f'],
        fitted_params_quinpirole_dict['tau_d'],
        fitted_params_quinpirole_dict['tau_rec'],
        fitted_params_quinpirole_dict['tau_Ca'],
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
    plt.savefig(f"/Users/reva/Documents/Python/SPSynapse/results/plots/trace/VL_{i}.pdf")


df_control = pd.DataFrame(fitted_params_control_list)
df_control['Condition'] = '10uM QP'
df_quinpirole = pd.DataFrame(fitted_params_quinpirole_list)
df_quinpirole['Condition'] = 'Quinpirole'

df_all = pd.concat([df_control, df_quinpirole], ignore_index=True)
df_all.to_csv("DF_All_VL.csv")
# List of parameters to analyze
parameters = ['calcium_channel_open_prob',  'total_Ca', 'total_P']

# Loop over each region and plot
for region in Regions:
    fig, axes = plt.subplots(len(parameters), 1, figsize=(10, len(parameters)*4))
    region_data = df_all[df_all['Region'] == region]

    for i, param in enumerate(parameters):
        sns.violinplot(data=region_data, x='Condition', y=param, ax=axes[i], inner=None, palette="Set2")
        sns.scatterplot(data=region_data, x='Condition', y=param, ax=axes[i], color="black", s=50)
        axes[i].set_title(f'{param} - {region}')

        # Statistical test (Mann-Whitney U test)
        control_data = region_data[region_data['Condition'] == '10uM QP'][param]
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
    plt.savefig(f"Stat_{region}.pdf")


# Create the figure
plt.figure(figsize=(8, 6))
for region in Regions:
    fig, axes = plt.subplots(len(parameters), 1, figsize=(10, len(parameters)*4))
    region_data = df_all[df_all['Region'] == region]

    for i, param in enumerate(parameters):
# Create half-violin plot for Baseline
        control_data = region_data[region_data['Condition'] == '10uM QP'][param]
        quinpirole_data = region_data[region_data['Condition'] == 'Quinpirole'][param]

        # Create half-violin plots
        sns.violinplot(data=region_data[region_data['Condition'] == '10uM QP'], x='Condition', y=param, ax=axes[i], 
                       inner=None, color='gray', alpha=0.3)
        sns.violinplot(data=region_data[region_data['Condition'] == 'Quinpirole'], x='Condition', y=param, ax=axes[i], 
                       inner=None, color='red', alpha=0.3)

        # Connect paired points
        for j in range(len(control_data)):
            axes[i].plot([0, 1], [control_data.iloc[j], quinpirole_data.iloc[j]], color='red', alpha=0.5)

        # Paired statistical test
        # Paired statistical test
        try:
            stat, p_value = wilcoxon(control_data, quinpirole_data, zero_method='zsplit')
        except ValueError as e:
            # HanVLe case when all pairs have zero difference
            if "x - y is zero for all elements" in str(e):
                p_value = 1.0  # Non-significant by default
            else:
                raise e

        if p_value < 0.001:
            significance = '***'
        elif p_value < 0.01:
            significance = '**'
        elif p_value < 0.05:
            significance = '*'
        else:
            significance = 'ns'  # not significant

        # Add statistical significance
        axes[i].text(0.5, 1.05, f'p = {p_value:.4f} {significance}', ha='center', va='bottom', 
                     transform=axes[i].transAxes, fontsize=12)

        axes[i].set_ylabel("Measurement Value")

    plt.tight_layout()
    plt.savefig(f"Stat_{region}_pairwise.pdf")
