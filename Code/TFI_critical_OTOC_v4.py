"""
TFI_critical_OTOC_v4.py — CONVENCIÓN EXACTA PROYECTO OTOC
- C(t) = |Tr(W(t) V W(t) V)| / dim  (infinite temperature)
- W = σ_z sitio 0, V = σ_z sitio 1 (r=1 FIJO, no r=N/2)
- Tiempos: t ∈ [0, D_MAX] con D_MAX = 5*N (convención proyecto)
- Ω = <|C(t)|>_late / C(0), late = segunda mitad de tiempos
- N = 4,5,6,7,8,9,10,11,12 (pares E impares)

TFI: H = -J Σ σ^z_i σ^z_{i+1} - h Σ σ^x_i, h/J=1 (crítico)
PBC (condiciones de frontera periódicas)

Nota: TFI puro es integrable. El científico de OTOC predice
anti-scrambling (Ω crece con N). Este script verifica.
También corre TFI + h_z (no integrable) para comparar.
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
import json, time, math

sx = np.array([[0,1],[1,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def kron_chain(ops):
    r = ops[0]
    for i in range(1,len(ops)): r = np.kron(r, ops[i])
    return r

def build_TFI(N, J=1.0, h_x=1.0, h_z=0.0):
    """
    H = -J Σ σ^z_i σ^z_{i+1} - h_x Σ σ^x_i - h_z Σ σ^z_i
    h_z = 0: integrable (Jordan-Wigner)
    h_z > 0: no integrable (scrambles)
    """
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(N):
        j = (i+1) % N
        ops = [I2]*N; ops[i] = sz; ops[j] = sz
        H -= J * kron_chain(ops)
    for i in range(N):
        ops = [I2]*N; ops[i] = sx
        H -= h_x * kron_chain(ops)
    if h_z != 0:
        for i in range(N):
            ops = [I2]*N; ops[i] = sz
            H -= h_z * kron_chain(ops)
    return H

def build_op(N, site):
    ops = [I2]*N; ops[site] = sz
    return kron_chain(ops)

def compute_OTOC_project(H, W, V, D_MAX, n_times=200):
    """
    Convención exacta del proyecto OTOC:
    C(t) = |Tr(W(t) V W(t) V)| / dim
    D_MAX = 5*N (máximo tiempo)
    """
    dim = H.shape[0]
    times = np.linspace(0.1, D_MAX, n_times)
    
    C0 = abs(np.trace(W @ V @ W @ V)) / dim
    
    C_t = np.zeros(n_times)
    for idx, t in enumerate(times):
        eiHt = expm(1j * H * t)
        Wt = eiHt @ W @ eiHt.conj().T
        C_t[idx] = abs(np.trace(Wt @ V @ Wt @ V)) / dim
    
    C_late = np.mean(C_t[n_times//2:])
    Omega = C_late / C0 if C0 > 1e-15 else 0
    
    return Omega, C0, C_late

def omega_fit(N, A, c):
    return 1.0 / (1.0 + A * N**c)

def run_model(name, N_values, J=1.0, h_x=1.0, h_z=0.0):
    print(f"\n{'='*60}")
    print(f"{name}")
    print(f"J={J}, h_x={h_x}, h_z={h_z}")
    print(f"r=1 (fijo), D_MAX=5N")
    print(f"{'='*60}")
    
    Omega_list = []
    N_used = []
    
    for N in N_values:
        dim = 2**N
        mem_gb = (dim**2 * 16) / (1024**3)
        
        if mem_gb > 6:
            print(f"N={N}: SKIP (mem={mem_gb:.2f}GB)")
            continue
        
        D_MAX = 5 * N  # Convención proyecto OTOC
        
        t0 = time.time()
        H = build_TFI(N, J, h_x, h_z)
        
        # r=1: W en sitio 0, V en sitio 1
        W = build_op(N, 0)
        V = build_op(N, 1)
        
        Omega, C0, C_late = compute_OTOC_project(H, W, V, D_MAX)
        dt = time.time() - t0
        
        print(f"N={N:2d}: Ω={Omega:.6f} C0={C0:.4f} C_late={C_late:.6f} ({dt:.0f}s)")
        
        Omega_list.append(Omega)
        N_used.append(N)
    
    # Análisis
    print(f"\nN:  {N_used}")
    print(f"Ω:  {[f'{o:.6f}' for o in Omega_list]}")
    
    # ¿Monótono decreciente?
    mono_dec = all(Omega_list[i] >= Omega_list[i+1] for i in range(len(Omega_list)-1))
    mono_inc = all(Omega_list[i] <= Omega_list[i+1] for i in range(len(Omega_list)-1))
    print(f"Monotono ↓: {mono_dec} (scrambling)")
    print(f"Monotono ↑: {mono_inc} (anti-scrambling)")
    
    # Fit si es decreciente
    result = {"name": name, "N": N_used, "Omega": Omega_list, 
              "h_z": h_z, "monotone_decreasing": mono_dec}
    
    if mono_dec and len(N_used) >= 3 and Omega_list[0] > Omega_list[-1] + 0.01:
        try:
            N_arr = np.array(N_used, dtype=float)
            Om_arr = np.array(Omega_list)
            po, pc = curve_fit(omega_fit, N_arr, Om_arr,
                              p0=[0.1, 1.0], bounds=([1e-10, 0.01], [1e6, 20]),
                              maxfev=50000)
            A, c = po
            c_err = np.sqrt(pc[1,1])
            Om_pred = omega_fit(N_arr, *po)
            SS_res = np.sum((Om_arr - Om_pred)**2)
            SS_tot = np.sum((Om_arr - np.mean(Om_arr))**2)
            R2 = 1 - SS_res/SS_tot if SS_tot > 1e-15 else 0
            
            print(f"\nFIT: c = {c:.4f} ± {c_err:.4f}, R² = {R2:.6f}")
            result["c"] = c
            result["c_err"] = c_err
            result["R2"] = R2
            
            # Comparaciones
            comps = [("ln√2",0.3466),("ln2",0.6931),("1/2",0.5),("1",1.0),
                     ("2",2.0),("π/2",math.pi/2),("π",math.pi),("√2",1.4142),
                     ("e",math.e),("φ",1.6180),("4/3",1.3333),("3/2",1.5),
                     ("16/7",2.2857),("2+√2",3.4142),("2π",6.2832),
                     ("3",3.0),("4",4.0),("5",5.0),("6",6.0)]
            print("Closest:")
            for n,v in sorted(comps, key=lambda x: abs(c-x[1]))[:5]:
                disc = abs(c-v)/v*100
                m = " ◄" if disc < 5 else ""
                print(f"  c={c:.4f} vs {n}={v:.4f}: {disc:.1f}%{m}")
        except Exception as e:
            print(f"Fit failed: {e}")
    elif not mono_dec:
        print("No fit — Ω no es monótono decreciente")
    
    return result

# ============================================================
# EJECUCIÓN PRINCIPAL
# ============================================================

all_results = []

# 1. TFI puro (integrable) — debería NO scrambliar
print("*" * 60)
print("PARTE 1: TFI PURO (INTEGRABLE)")
print("Predicción del científico OTOC: anti-scrambling")
print("*" * 60)
r1 = run_model("TFI puro h/J=1 (integrable)", 
               [4, 5, 6, 7, 8, 9, 10, 11, 12],
               J=1.0, h_x=1.0, h_z=0.0)
all_results.append(r1)

# 2. TFI + h_z = 0.1 (débilmente no integrable)
print("\n" + "*" * 60)
print("PARTE 2: TFI + h_z (NO INTEGRABLE)")
print("Predicción: scrambling emerge con h_z")
print("*" * 60)

for h_z in [0.1, 0.3, 0.5]:
    r = run_model(f"TFI + h_z={h_z}",
                  [4, 5, 6, 7, 8, 9, 10],  # hasta 10 para velocidad
                  J=1.0, h_x=1.0, h_z=h_z)
    all_results.append(r)

# Guardar
with open("TFI_v4_all_results.json", 'w') as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\n{'='*60}")
print("Saved TFI_v4_all_results.json")
print(f"{'='*60}")
