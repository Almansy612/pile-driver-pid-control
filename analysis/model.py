"""Pile-driver rig: active vibration damping with a hydraulic actuator + PID.

Single source of truth for every number in the README. Run:
    python analysis/model.py            # prints parameters, tuning, results
    python analysis/model.py --plots    # also regenerates images/*.png
"""
from __future__ import annotations
import math, sys
import numpy as np

# ── 1. Rig / pile structure (single-degree-of-freedom) ──────────────────────
M      = 1000.0                 # kg   effective rig mass
K      = 2.0e6                  # N/m  effective stiffness
ZETA   = 0.05                   # –    structural damping ratio (steel rig, typical 2–5 %)
WN     = math.sqrt(K / M)       # rad/s
FN     = WN / (2 * math.pi)     # Hz
C      = 2 * ZETA * math.sqrt(K * M)   # N·s/m

# ── 2. Hydraulic actuator (double-acting cylinder, valve-controlled) ────────
D_BORE = 0.050                  # m
D_ROD  = 0.020                  # m
STROKE = 0.200                  # m
P_MAX  = 200e5                  # Pa  (200 bar)
BETA   = 1.5e9                  # Pa  oil bulk modulus (≈1.5 GPa for mineral oil)
Q_MAX  = 20 / 60_000            # m³/s (20 L/min)
A_BORE = math.pi / 4 * D_BORE**2
A_ANN  = math.pi / 4 * (D_BORE**2 - D_ROD**2)   # rod side: the weaker direction
F_MAX  = P_MAX * A_ANN          # N   symmetric force limit used by the controller
V_MID  = A_ANN * STROKE / 2     # m³  chamber volume at mid-stroke (linearisation point)
# Chamber pressure: dP/dt = (β/V)(Kq·u − Kc·P), with Kq = Q_max and Kc = Q_max/P_max
# so that u = 1 settles at P_max. That is first order with
TAU_A  = V_MID * P_MAX / (BETA * Q_MAX)          # s
# → F(s)/u(s) = F_MAX / (TAU_A s + 1),  u ∈ [−1, 1] normalised valve command

# ── 3. Disturbance: ram impact reaction reaching the rig ────────────────────
F_IMPACT = 20e3                 # N   half-sine peak
T_IMPACT = 0.025                # s   pulse duration
J_IMPACT = 2 / math.pi * F_IMPACT * T_IMPACT     # N·s impulse

def impact(t):
    return F_IMPACT * math.sin(math.pi * t / T_IMPACT) if 0.0 <= t < T_IMPACT else 0.0

# ── 4. Controller: PID on rig displacement (setpoint 0), filtered derivative ─
#      u = −(Kp·x + Ki·∫x dt + Kd·ẋ_f),  ẋ_f = N·s/(s+N) · x,  then clipped to ±1
TUNED    = dict(Kp=20.0, Ki=1200.0, Kd=1.8, N=600.0)
#   From the grid search in analysis/tune.py (settling across the impact AND
#   static-load tests, subject to PM ≥ 45°, GM ≥ 10 dB, |u| ≤ 0.8). Within 4 ms of
#   the fastest candidate; preferred for its larger gain margin (21 vs 19 dB) and
#   lower impact peak.
W_RAM = 1000 * 9.81             # N  static load step: ram weight lifted onto the rig
ORIGINAL = {  # the three gain sets from the coursework report
    "Report 'ideal case'": dict(Kp=50_000, Ki=1_000, Kd=5_000, N=100.0),
    "Report trial 1":      dict(Kp=199_000, Ki=10_000, Kd=1_000, N=100.0),
    "Report trial 2":      dict(Kp=300_000, Ki=50_000, Kd=0.0,  N=50.0),
}

# ── Simulation (nonlinear: includes actuator saturation) ────────────────────
def simulate(gains=None, t_end=2.0, dt=1e-5, saturate=True, dist=None):
    """States: x, v, F (actuator force), xi (∫x), xf (derivative filter state)."""
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    X = np.zeros((n, 3)); U = np.zeros(n)
    dist = dist or impact
    x = v = F = xi = xf = 0.0
    for i in range(n):
        X[i] = (x, v, F)
        if gains:
            N = gains["N"]
            xdot_f = N * (x - xf)                     # filtered derivative of x
            u = -(gains["Kp"] * x + gains["Ki"] * xi + gains["Kd"] * xdot_f)
            if saturate: u = max(-1.0, min(1.0, u))
        else:
            u = 0.0; N = 1.0; xdot_f = 0.0
        U[i] = u
        # semi-implicit Euler at 10 µs (actuator τ ≈ 6.6 ms, rig period ≈ 140 ms)
        a = (dist(t[i]) + F - C * v - K * x) / M
        v += a * dt; x += v * dt
        F += (F_MAX * u - F) / TAU_A * dt
        xi += x * dt
        if gains: xf += N * (x - xf) * dt
    return t, X[:, 0], X[:, 1], X[:, 2], U

def settle_time(t, y, band):
    out = np.where(np.abs(y) > band)[0]
    return float(t[out[-1]]) if len(out) else 0.0

def static_load(t):
    return W_RAM

def metrics(gains=None, t_end=2.0, dist=None):
    t, x, v, F, u = simulate(gains, t_end, dist=dist)
    peak = float(np.max(np.abs(x)))
    return dict(peak_mm=peak * 1e3,
                settle_s=settle_time(t, x, 0.1e-3),          # until |x| < 0.1 mm for good
                rms_mm=float(np.sqrt(np.mean(x**2))) * 1e3,
                u_peak=float(np.max(np.abs(u))),
                F_peak_kN=float(np.max(np.abs(F))) / 1e3,
                final_mm=float(x[-1]) * 1e3)

# ── Linear loop analysis (python-control) ───────────────────────────────────
def loop_tf(g):
    import control as ct
    s = ct.tf("s")
    Cpid = g["Kp"] + g["Ki"] / s + g["Kd"] * g["N"] * s / (s + g["N"])
    Act = F_MAX / (TAU_A * s + 1)
    Plant = 1 / (M * s**2 + C * s + K)
    return Cpid * Act * Plant, Plant

def margins(g):
    import control as ct
    L, _ = loop_tf(g)
    gm, pm, wcg, wcp = ct.margin(L)
    cl = ct.feedback(L, 1)
    poles = ct.poles(cl)
    return dict(gm_db=20 * math.log10(gm) if gm and np.isfinite(gm) else float("inf"),
                pm_deg=pm, wc_rad_s=wcp, stable=bool(np.all(poles.real < 0)),
                zeta_min=float(min(-p.real / abs(p) for p in poles if abs(p) > 1e-9)))

def main():
    global TAU_A, M
    print("== Parameters")
    print(f"rig: m={M:.0f} kg  k={K:.2e} N/m  ζ={ZETA}  → ωn={WN:.2f} rad/s  fn={FN:.2f} Hz  c={C:.0f} N·s/m")
    print(f"cylinder: A_bore={A_BORE*1e3:.3f}e-3 m²  A_annulus={A_ANN*1e3:.3f}e-3 m²  F_max={F_MAX/1e3:.2f} kN")
    print(f"V_mid={V_MID*1e6:.1f} cm³  τ_a = V·P/(β·Q) = {TAU_A*1e3:.2f} ms  (static gain {F_MAX/1e3:.2f} kN per unit command)")
    print(f"impact: half-sine {F_IMPACT/1e3:.0f} kN × {T_IMPACT*1e3:.0f} ms  (J = {J_IMPACT:.0f} N·s)")
    print("\n== Loop margins, tuned PID", TUNED)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in margins(TUNED).items()})
    print("\n== Impact response (nonlinear, with ±F_max saturation)")
    rows = [("Uncontrolled", None), ("Tuned PID", TUNED)] + list(ORIGINAL.items())
    print(f"{'case':24s} {'peak mm':>8s} {'settle s':>9s} {'rms mm':>7s} {'|u|max':>7s} {'F kN':>6s}  stable")
    for name, g in rows:
        m = metrics(g)
        st = "—" if g is None else margins(g)["stable"]
        print(f"{name:24s} {m['peak_mm']:8.2f} {m['settle_s']:9.3f} {m['rms_mm']:7.3f} {m['u_peak']:7.2f} {m['F_peak_kN']:6.1f}  {st}")
    print("\n== Static load: ram weight (9.81 kN) applied as a step")
    for name, g in [("Uncontrolled", None), ("Tuned PID", TUNED)]:
        mm = metrics(g, dist=static_load)
        print(f"{name:24s} peak {mm['peak_mm']:.2f} mm   final {mm['final_mm']:.3f} mm   |u|max {mm['u_peak']:.2f}")
    print("\n== Robustness of the tuned loop (linear margins)")
    base = (TAU_A, M)
    for label, tau, mass in [("nominal", base[0], base[1]), ("τ_a ×2 (full-stroke volume)", base[0]*2, base[1]),
                             ("mass −20 %", base[0], base[1]*0.8), ("mass +20 %", base[0], base[1]*1.2)]:
        TAU_A, M = tau, mass
        mg = margins(TUNED)
        print(f"{label:30s} PM {mg['pm_deg']:5.1f}°  GM {mg['gm_db']:5.1f} dB  stable {mg['stable']}")
    TAU_A, M = base
    if "--plots" in sys.argv: plots()

def plots():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import control as ct
    RED, BLUE, GREY = "#c0392b", "#1f6fb2", "#7f8c8d"

    for fname, dist, title in [
        ("impact-response.png", impact, "Ram-impact pulse (20 kN half-sine, 25 ms) hitting the rig"),
        ("static-load.png", static_load, "Ram weight (9.81 kN) transferred onto the rig as a step")]:
        t, x0, *_ = simulate(None, dist=dist)
        t, x1, v1, F1, u1 = simulate(TUNED, dist=dist)
        fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 2]))
        ax[0].plot(t, x0 * 1e3, color=RED, lw=1.2, label="Uncontrolled (ζ = 0.05)")
        ax[0].plot(t, x1 * 1e3, color=BLUE, lw=2.0, label="Hydraulic actuator + PID")
        ax[0].axhline(0, color=GREY, lw=.8)
        ax[0].set_ylabel("rig displacement (mm)"); ax[0].legend(loc="upper right"); ax[0].grid(alpha=.3)
        ax[0].set_title(title)
        ax[1].plot(t, F1 / 1e3, color=BLUE, lw=1.6)
        for sgn in (1, -1): ax[1].axhline(sgn * F_MAX / 1e3, ls="--", color=GREY, lw=1)
        ax[1].text(t[-1], F_MAX / 1e3 * 0.8, "actuator limit ±33 kN (200 bar)", ha="right", color=GREY, fontsize=9)
        ax[1].set_ylabel("actuator force (kN)"); ax[1].set_xlabel("time (s)"); ax[1].grid(alpha=.3)
        ax[1].set_xlim(0, t[-1])
        fig.tight_layout(); fig.savefig("images/" + fname, dpi=140); plt.close(fig)

    L, P = loop_tf(TUNED)
    w = np.logspace(0, 4, 1000); f = w / (2 * math.pi)
    rP = ct.frequency_response(P, w); rD = ct.frequency_response(ct.feedback(P, L / P), w)
    rL = ct.frequency_response(L, w)
    mg = margins(TUNED)
    fig, ax = plt.subplots(3, 1, figsize=(9, 8.5), sharex=True)
    ax[0].semilogx(f, 20 * np.log10(np.abs(rP.complex) * K), color=RED, label="Uncontrolled")
    ax[0].semilogx(f, 20 * np.log10(np.abs(rD.complex) * K), color=BLUE, lw=2, label="With PID")
    ax[0].set_ylabel("|x / F_d| · k  (dB)"); ax[0].legend(loc="lower left"); ax[0].grid(alpha=.3, which="both")
    ax[0].set_title("Disturbance force → rig displacement: 7.1 Hz resonance flattened")
    ax[1].semilogx(f, 20 * np.log10(np.abs(rL.complex)), color=BLUE); ax[1].axhline(0, color=GREY, lw=1)
    ax[1].set_ylabel("|L(jω)| (dB)"); ax[1].grid(alpha=.3, which="both")
    ax[1].set_title(f"Open loop L = PID · actuator · rig   —   PM {mg['pm_deg']:.0f}°,  GM {mg['gm_db']:.0f} dB,  crossover {mg['wc_rad_s']/(2*math.pi):.0f} Hz")
    ax[2].semilogx(f, np.degrees(np.unwrap(np.angle(rL.complex))), color=BLUE); ax[2].axhline(-180, color=GREY, lw=1)
    ax[2].set_ylabel("∠L(jω) (deg)"); ax[2].set_xlabel("frequency (Hz)"); ax[2].grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig("images/frequency-response.png", dpi=140); plt.close(fig)

    # original coursework gains on the corrected plant
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for (name, g), col in zip(ORIGINAL.items(), ["#8e44ad", "#d35400", "#16a085"]):
        t, x, *_ = simulate(g, t_end=1.0)
        ax.plot(t, x * 1e3, color=col, lw=1.2, label=name)
    t, x, *_ = simulate(TUNED, t_end=1.0)
    ax.plot(t, x * 1e3, color=BLUE, lw=2, label="Re-tuned PID")
    ax.set_ylim(-60, 60); ax.set_xlabel("time (s)"); ax.set_ylabel("displacement (mm)"); ax.grid(alpha=.3)
    ax.set_title("Original report gains on the corrected model: all three go unstable and saturate")
    ax.legend(loc="lower left", fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig("images/original-gains.png", dpi=140); plt.close(fig)
    print("plots written to images/")

if __name__ == "__main__":
    main()
