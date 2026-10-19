"""Grid search used to pick the PID gains in model.py.

Objective: shortest settling time across BOTH tests (impact pulse to |x| < 0.1 mm,
ram-weight step to |x| < 0.1 mm), subject to
    phase margin ≥ 45°, gain margin ≥ 10 dB, peak valve command |u| ≤ 0.8.
Linear models (python-control) for speed; model.py re-checks the winner with
the nonlinear, saturating simulation.
"""
import itertools, math, numpy as np, control as ct
import model as m

s = ct.tf("s")
t = np.linspace(0, 2, 6001)
Fd = np.array([m.impact(x) for x in t])
Fs = np.full_like(t, m.W_RAM)
P = 1 / (m.M * s**2 + m.C * s + m.K)
A = m.F_MAX / (m.TAU_A * s + 1)

rows = []
for Kp, Kd, Ki, N in itertools.product([0, 10, 20, 40, 60, 90, 120],
                                       [1.0, 1.4, 1.8, 2.2, 2.6, 3.0],
                                       [0, 200, 400, 800, 1200],
                                       [300, 600]):
    C = Kp + Ki / s + Kd * N * s / (s + N)
    gm, pm, _, _ = ct.margin(C * A * P)
    gm_db = 20 * math.log10(gm) if np.isfinite(gm) else float("inf")
    if pm < 45 or gm_db < 10:
        continue
    T = ct.feedback(P, C * A)
    _, x = ct.forced_response(T, t, Fd)
    _, u = ct.forced_response(-C * T, t, Fd)
    _, xs = ct.forced_response(T, t, Fs)
    if np.max(np.abs(u)) > 0.8:
        continue
    ts = t[np.where(np.abs(x) > 1e-4)[0][-1]]
    tss = t[np.where(np.abs(xs) > 1e-4)[0][-1]]
    rows.append((max(ts, tss), ts, tss, pm, gm_db, np.max(np.abs(u)), Kp, Kd, Ki, N))

rows.sort()
print(" worst_ts  impact_ts static_ts    PM     GM   |u|   Kp   Kd    Ki    N")
for r in rows[:10]:
    print(f"{r[0]:8.3f} {r[1]:9.3f} {r[2]:9.3f} {r[3]:6.1f} {r[4]:6.1f} {r[5]:5.2f} {r[6]:4} {r[7]:4} {r[8]:5} {r[9]:4}")
