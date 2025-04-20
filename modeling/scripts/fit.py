import numpy as np
from model import synaptic_response
from config import PPR_data

control_params_dict = {}  # will be filled per cell

# Objective function for control condition
def objective_function(params, condition='control'):
    observed_PPR = PPR_data[condition]
    PPR_model, IPSC_model, _, _, _, _ = synaptic_response(params, condition)
    error = (PPR_model[0] - observed_PPR) ** 2
    return error

# Objective function for quinpirole with full parameter optimization
def objective_function_quinpirole(params):
    P0, calcium_influx, tau_f, tau_d, tau_Ca, calcium_channel_open_prob = params

    Kc = control_params_dict['Kc']
    tau_rec = control_params_dict['tau_rec']

    model_params_q = [P0, Kc, tau_f, tau_d, tau_rec, tau_Ca, calcium_influx, calcium_channel_open_prob]
    model_params_c = [
        control_params_dict['P0'],
        Kc,
        control_params_dict['tau_f'],
        control_params_dict['tau_d'],
        tau_rec,
        control_params_dict['tau_Ca'],
        control_params_dict['calcium_influx'],
        control_params_dict['calcium_channel_open_prob']
    ]

    PPR_model_q, IPSC_q, _, _, _, _ = synaptic_response(model_params_q, 'Quinpirole')
    PPR_model_c, IPSC_c, _, _, _, _ = synaptic_response(model_params_c, 'control')

    peak_q = np.max(IPSC_q)
    peak_c = np.max(IPSC_c)
    observed_PPR_q = PPR_data['Quinpirole']

    error = (PPR_model_q[0] - observed_PPR_q) ** 2
    penalty = 0

    if peak_q > peak_c:
        penalty += (peak_q - peak_c)**2 * 1e3

    if peak_q != 0:
        ratio = peak_c / peak_q
    else:
        ratio = np.inf

    if ratio > 4:
        penalty += ratio ** 5 * 1e3

    return error + penalty

# Objective with scaled calcium channel open probability
def objective_function_quinpirole_p(scaling_params):
    observed_PPR_q = PPR_data['Quinpirole']

    model_params_c = [
        control_params_dict['P0'],
        control_params_dict['Kc'],
        control_params_dict['tau_f'],
        control_params_dict['tau_d'],
        control_params_dict['tau_rec'],
        control_params_dict['tau_Ca'],
        control_params_dict['calcium_influx'],
        control_params_dict['calcium_channel_open_prob']
    ]

    PPR_model_q, IPSC_q, _, _, _, _ = synaptic_response(
        model_params_c, 'Quinpirole', scaling_params=scaling_params
    )

    error = (PPR_model_q[0] - observed_PPR_q) ** 2

    _, IPSC_c, _, _, _, _ = synaptic_response(model_params_c, 'control')
    peak_q = np.max(IPSC_q)
    peak_c = np.max(IPSC_c)

    ratio = peak_c / peak_q if peak_q != 0 else np.inf
    penalty = (ratio ** 5 * 1e3) if ratio > 4 else 0

    return error + penalty