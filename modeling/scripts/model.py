import numpy as np
from config import time, dt, spike_times, Ca0, m, n_infty, tau_ipsc

def synaptic_response(params, condition, scaling_params=None):
    P0, Kc, tau_f, tau_d, tau_rec, tau_Ca, calcium_influx, calcium_channel_open_prob = params

    u = np.zeros_like(time)
    R = np.ones_like(time)
    Ca = np.zeros_like(time)
    P = np.zeros_like(time)
    A = np.zeros_like(time)
    n = np.ones_like(time)
    IPSC = np.zeros_like(time)
    Ca[0] = Ca0

    u[0] = P0
    R[0] = 1
    n[0] = n_infty

    if scaling_params is not None and condition == 'Quinpirole':
        scale_open_prob, scale_d, scale_f = scaling_params
        calcium_channel_open_prob *= scale_open_prob

    for t in range(1, len(time)):
        dCa_dt = -(Ca[t-1] - Ca0) / tau_Ca
        if any(np.isclose(time[t], spike_times, atol=dt/2)):
            dCa_dt += (calcium_influx * calcium_channel_open_prob) / tau_Ca

        Ca[t] = Ca[t-1] + dCa_dt * dt
        P[t] = P0 * (Ca[t]**m / (Kc**m + Ca[t]**m))
        P[t] = np.clip(P[t], 0, 1)

        dn = (n_infty - n[t-1]) / tau_rec * dt
        if time[t] in spike_times:
            dn -= u[t-1] * R[t-1] * n[t-1]
        n[t] = np.clip(n[t-1] + dn, 0, n_infty)

        du = -u[t-1] / tau_f * dt
        if time[t] in spike_times:
            du += P[t] * (1 - u[t-1])
        u[t] = np.clip(u[t-1] + du, 0, 1)

        dR = (1 - R[t-1]) / tau_d * dt
        if time[t] in spike_times:
            dR -= u[t-1] * R[t-1] * n[t-1]
        R[t] = np.clip(R[t-1] + dR, 0, 1)

        if time[t] in spike_times:
            A[t] = u[t] * R[t] * n[t]
            IPSC[t:] += A[t] * np.exp(-(time[t:] - time[t]) / tau_ipsc)

    idx1 = np.where(np.isclose(time, spike_times[0], atol=dt/2))[0]
    idx2 = np.where(np.isclose(time, spike_times[1], atol=dt/2))[0]

    if idx1.size > 0 and idx2.size > 0:
        A1 = A[idx1[0]]
        A2 = A[idx2[0]]
        PPR_model = A2 / A1 if A1 != 0 else np.inf
    else:
        PPR_model = np.nan

    total_P = np.trapz(P, time)
    total_Ca = np.trapz(Ca, time)

    return np.array([PPR_model]), IPSC, P, Ca, total_P, total_Ca