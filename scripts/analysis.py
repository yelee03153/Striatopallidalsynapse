import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from model import synaptic_response
from fit import objective_function, objective_function_quinpirole_p, control_params_dict
from config import PPR_data, param_bounds

# Load experimental paired-pulse and amplitude data
all_data = pd.read_pickle("results/features/features_all_data.pkl")


def run_fitting_pipeline(region, n_fits_required=10, cond="control"):
    """
    Run model fitting for a specific brain region.
    Returns a DataFrame with fitted control and quinpirole parameters.
    """
    #Type of experiment
    if cond=="control":
        ex_type="10uM QP"
    elif cond=="6OHDA":
        ex_type="6-OHDA"

    # Extract relevant experimental data for the selected region
    ppr_values = all_data["PPR"].loc[ex_type].loc[region]["Baseline"].dropna()
    ppr_values_quinpirole = all_data["PPR"].loc[ex_type].loc[region]["During"].dropna()

    # Store results
    fitted_params_control_list = []
    fitted_params_quinpirole_list = []
    attempts = 0
    valid_fits = 0
    max_attempts = 100*n_fits_required

    while valid_fits < n_fits_required and  attempts < max_attempts:

        PPR_data['control'] = np.mean(ppr_values)
        PPR_data['Quinpirole'] = np.mean(ppr_values_quinpirole)

        result_control = differential_evolution(objective_function, param_bounds, strategy='best1bin', maxiter=50, seed=np.random.randint(0, 1e6))
        control_params = result_control.x

        _, _, _, _, tot_P_c, tot_Ca_c = synaptic_response(control_params, 'control')

        control_params_dict.update({
            'Region': region,
            'P0': control_params[0],
            'Kc': control_params[1],
            'tau_f': control_params[2],
            'tau_d': control_params[3],
            'tau_rec': control_params[4],
            'tau_Ca': control_params[5],
            'calcium_influx': control_params[6],
            'calcium_channel_open_prob': control_params[7],
            'total_P': tot_P_c,
            'total_Ca': tot_Ca_c
        })

        bounds_scaling = [(0.01, 1.01), (0.8, 1.2), (0.8, 1.2)]
        result_quinpirole = differential_evolution(objective_function_quinpirole_p, bounds_scaling, strategy='best1bin', maxiter=150, seed=np.random.randint(0, 1e6))
        scaling_params_quinpirole = result_quinpirole.x

        PPR_model_quinpirole, IPSC_quinpirole, _, _, tot_P_q, tot_Ca_q = synaptic_response(
            control_params, 'Quinpirole', scaling_params=scaling_params_quinpirole
        )

        if (PPR_model_quinpirole[0] > 1 and PPR_data['Quinpirole'] < 1) or \
        (PPR_model_quinpirole[0] < 1 and PPR_data['Quinpirole'] > 1) or \
        abs(PPR_model_quinpirole[0] - PPR_data['Quinpirole']) > 0.4:
            attempts += 1
            continue  # Skip this fit

        fitted_params_control_list.append(control_params_dict.copy())
        fitted_params_quinpirole_list.append({
            'Region': region,
            'P0': control_params[0],
            "Kc": control_params[1],
            'tau_f': control_params[2],
            'tau_d': control_params[3],
            'tau_rec': control_params[4],
            'tau_Ca': control_params[5],
            'calcium_influx': control_params[6],
            'calcium_channel_open_prob': control_params[7] * scaling_params_quinpirole[0],
            'total_P': tot_P_q,
            'total_Ca': tot_Ca_q
        })

        valid_fits += 1
        attempts += 1


    # Combine into DataFrames
    fitted_df_control = pd.DataFrame(fitted_params_control_list)
    fitted_df_control['Condition'] = 'Baseline'

    fitted_df_quinpirole = pd.DataFrame(fitted_params_quinpirole_list)
    fitted_df_quinpirole['Condition'] = 'Quinpirole'

    df_all = pd.concat([fitted_df_control, fitted_df_quinpirole], ignore_index=True)

    df_summary = df_all[['Region', 'Condition', 'calcium_channel_open_prob', 'total_Ca']]

    return df_summary

