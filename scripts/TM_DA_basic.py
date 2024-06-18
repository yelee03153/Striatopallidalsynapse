import numpy as np
import matplotlib.pyplot as plt

# Parameters
Ca0 = 0.01
Kd_DA = 10e-9  # Dissociation constant for D2 receptor (10 nM converted to M)
n_Hill = 1  # Hill coefficient for DA effect on calcium
P0 = 0.1
Kc = 0.5  # Dissociation constant for calcium effect on release probability
m = 4  # Hill coefficient for calcium effect on release probability
tau_f = 0.05
tau_d = 0.02
tau_DA = 0.05
tau_rec = 0.05  # Time constant for vesicle replenishment
dt = 0.001
time = np.arange(0, 10, dt)
DA_concentration = 0.5  # Stable extracellular DA concentration (can be set as needed)
n_infty = 1.0  # Steady-state vesicle availability
spike_times = [1, 1.2, 5, 7, 9]

# Variables
u = np.zeros_like(time)
R = np.ones_like(time)
Ca = np.zeros_like(time)
P = np.zeros_like(time)
A = np.zeros_like(time)
n = np.ones_like(time)  # Vesicle availability

# Initial conditions
u[0] = P0
R[0] = 1
n[0] = n_infty

for t in range(1, len(time)):
    # Nonlinear effect of DA on intracellular calcium concentration
    Ca[t] = Ca0 + (DA_concentration**n_Hill / (Kd_DA**n_Hill + DA_concentration**n_Hill))
    
    # Calcium-dependent release probability
    P[t] = P0 * (Ca[t]**m / (Kc**m + Ca[t]**m))
    
    # Update vesicle availability
    dn = (n_infty - n[t-1]) / tau_rec * dt
    if time[t] in spike_times:
        dn -= u[t-1] * R[t-1] * n[t-1]
    n[t] = n[t-1] + dn
    
    # Update u with facilitation
    du = (P[t] - u[t-1]) / tau_f * dt
    if time[t] in spike_times:
        du += P[t] * (1 - u[t-1])
    u[t] = u[t-1] + du
    
    # Update R with depression
    dR = (1 - R[t-1]) / tau_d * dt
    if time[t] in spike_times:
        dR -= u[t-1] * R[t-1] * n[t-1]
    R[t] = R[t-1] + dR
    
    # Synaptic response
    if time[t] in spike_times:
        A[t] = u[t] * R[t] * n[t]

# Plot results
plt.figure(figsize=(12, 10))
plt.subplot(511)
plt.plot(time, u, label='Facilitation (u)')
plt.legend()
plt.subplot(512)
plt.plot(time, A, label='Synaptic Response (A)')
plt.legend()
plt.subplot(513)
plt.plot(time, P, label='Release Probability (P)')
plt.legend()
plt.subplot(514)
plt.plot(time, n, label='Vesicle Availability (n)')
plt.legend()
plt.subplot(515)
plt.plot(time, R, label='Depression (R)')
plt.legend()
plt.show()
