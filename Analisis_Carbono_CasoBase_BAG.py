#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
carbono_caso_base.py
====================
Modelo desacoplado de formacion y gestion de carbono para el reactor de
estano liquido (DRM). Equivalente en Python de la hoja de calculo Excel del
caso base (Seccion 4.2.3 de la tesis de B. Agun, UPM).

Recibe las presiones parciales que DWSIM (DRM + RWGS, cinetica L-H) produce
A CADA TEMPERATURA y devuelve, por temperatura:
    - velocidades de craqueo de CH4 y Boudouard,
    - velocidad neta de formacion de carbono,
    - cobertura superficial de carbono en estado estacionario (theta_C^SS),
    - tiempo hasta saturacion del metal (t_sat),
    - perfil temporal de acumulacion C_bulk(t) y del factor CF(t).

Soporta barrido en temperatura (p. ej. 723-1273 K cada 50 K). Las presiones
de cada T deben venir de una corrida de DWSIM a esa misma T.

Independiente de DWSIM: las presiones se introducen a mano o se leen de CSV.
Admite presiones [Pa], fracciones molares y_i o caudales molares F_i [mol/s].
Dependencias: estandar (math, csv). numpy y matplotlib opcionales.
"""

import math
import csv

try:
    import numpy as np
    _HAS_NP = True
except ImportError:
    _HAS_NP = False

try:
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    import pandas as pd        # solo para leer los .xlsx de DWSIM
    _HAS_PD = True
except ImportError:
    _HAS_PD = False


# =============================================================================
# 1. CONFIGURACION  ---  VERIFICAR CONTRA EL EXCEL DEL CASO BASE
# =============================================================================
R = 8.314  # J/(mol.K)

# ---- Geometria del reactor (Plevan et al. 2015) -----------------------------
D_INT   = 35.9e-3      # m    diametro interior
H_LM    = 600e-3       # m    altura de llenado de Sn
V_R     = 6.07e-4      # m3   volumen reactivo
RHO_CAT = 6195.0       # kg/m3  catalyst loading (drift-flux, Catalan & Rezaei)
RHO_SN  = 6400.0       # kg/m3  densidad del Sn liquido puro
W_SN    = 3.76         # kg

# Densidad usada en t_sat y acumulacion. La tesis escribe rho_Sn (=6400).
# AVISO: W_Sn/V_R = 3.76/6.07e-4 = 6195 = RHO_CAT. Si la velocidad esta
# referida a kg de catalyst loading, lo coherente es RHO_CAT (~3% en t_sat).
RHO_ACUM = RHO_SN

T_BASE  = 1173.15      # K    (900 C)
P_BASE  = 101325.0     # Pa

# ---- Escenarios de incertidumbre (k_diss, C_sat) ----------------------------
# k_diss solo interviene en theta_C^SS; las curvas C_bulk(t) difieren entre
# escenarios solo por C_sat. Emparejamiento por columna (CONFIRMAR).
SCENARIOS = {
    "conservador": dict(k_diss=1e-7, C_sat=5.0),
    "central":     dict(k_diss=1e-4, C_sat=20.0),   # caso base (Okamoto 2012)
    "optimista":   dict(k_diss=3e-2, C_sat=50.0),
}

SPECIES = ("CH4", "CO2", "CO", "H2", "H2O")


# =============================================================================
# 2. CONSTANTES CINETICAS  (dependen de T)
# =============================================================================
def kinetic_constants(T):
    """Constantes cineticas y de adsorcion a temperatura T [K].
    Adsorcion en Pa^-1 (bar^-1 / 1e5), identicas a DRM/RWGS de DWSIM.
    Velocidad en mol/(kg.s). Keq_crack en bar, Keq_Boud en bar^-1."""
    RT = R * T
    return {
        "K1":   2.60e-2 * math.exp( 24410.0 / RT) / 1e5,  # K_CH4
        "K2":   2.61e-2 * math.exp( 22585.0 / RT) / 1e5,  # K_CO2
        "K3":   1.00e-5 * math.exp( 30000.0 / RT) / 1e5,  # K_H2
        "Km21": 1.00e-5 * math.exp(  8220.0 / RT) / 1e5,  # K_CO   (K_{-21})
        "Km23": 4.73e-6 * math.exp( 58662.0 / RT) / 1e5,  # K_H2O  (K_{-23})
        "k_crack": 6.95e3  * math.exp(-153800.0 / RT),
        "k_Boud":  1.34e15 * math.exp(-445000.0 / RT),
        "Keq_crack": 2.98e5 * math.exp(-84400.0 / RT),           # bar
        "Keq_Boud":  1.0 / (1.93e9 * math.exp(-168527.0 / RT)),  # bar^-1
    }


# =============================================================================
# 3. DENOMINADOR OMEGA Y VELOCIDADES (locales, en un punto)
# =============================================================================
def omega(p, K):
    """Denominador L-H comun (identico a DRM/RWGS). Presiones en Pa."""
    return (1.0
            + K["K1"]   * p["CH4"]
            + K["K2"]   * p["CO2"]
            + K["Km21"] * p["CO"]
            + math.sqrt(K["K3"] * p["H2"])
            + K["Km23"] * p["H2O"])


def carbon_rates(p, T, K=None):
    """Velocidades LOCALES de carbono [mol/(kg.s)] en un punto.
    p: presiones parciales [Pa].  Prefactores L-H en Pa; corchetes de
    equilibrio en bar (Keq en bar / bar^-1)."""
    if K is None:
        K = kinetic_constants(T)
    Om = omega(p, K)
    Om2 = Om * Om
    b = {sp: p[sp] / 1e5 for sp in p}  # bar, solo para los corchetes

    if b["CH4"] > 0:                                   # CH4 -> C + 2 H2
        aff_crack = 1.0 - (b["H2"] ** 2) / (K["Keq_crack"] * b["CH4"])
    else:
        aff_crack = 0.0
    r_crack = K["k_crack"] * K["K1"] * p["CH4"] * aff_crack / Om2

    if b["CO"] > 0:                                    # 2 CO -> C + CO2
        aff_boud = 1.0 - b["CO2"] / (K["Keq_Boud"] * b["CO"] ** 2)
    else:
        aff_boud = 0.0
    r_Boud = K["k_Boud"] * K["Km21"] * p["CO"] * aff_boud / Om2

    return {"r_crack": r_crack, "r_Boud": r_Boud, "r_net": r_crack + r_Boud,
            "omega": Om, "aff_crack": aff_crack, "aff_boud": aff_boud}


# =============================================================================
# 4. PROMEDIO SOBRE EL PERFIL AXIAL
# =============================================================================
def mean_pressures(profile):
    """Media aritmetica de un perfil axial de presiones [Pa]."""
    keys = profile[0].keys()
    n = len(profile)
    return {k: sum(row[k] for row in profile) / n for k in keys}


def mean_carbon_rates(profile, T, K=None, mode="rates"):
    """Velocidades medias de carbono sobre el lecho [mol/(kg.s)].

    profile : lista de dicts de presiones [Pa] (puntos axiales). Tambien
              admite un unico dict (un punto, p. ej. la salida).
    mode = 'rates'     -> media de r(p_i) punto a punto  (RIGUROSO; r es no
                          lineal en p; recomendado).
         = 'pressures' -> r(media de p_i)  (equivalente al Excel de
                          'presiones medias').
    Nota: la media es aritmetica (asume puntos axiales de igual volumen).
    """
    if isinstance(profile, dict):
        profile = [profile]
    if K is None:
        K = kinetic_constants(T)
    if mode == "pressures":
        return carbon_rates(mean_pressures(profile), T, K)
    n = len(profile)
    acc = {"r_crack": 0.0, "r_Boud": 0.0, "r_net": 0.0,
           "omega": 0.0, "aff_crack": 0.0, "aff_boud": 0.0}
    for p in profile:
        r = carbon_rates(p, T, K)
        for k in acc:
            acc[k] += r[k] / n
    return acc


# =============================================================================
# 5. COBERTURA, TIEMPO DE SATURACION E INTEGRACION TEMPORAL
# =============================================================================
def theta_ss(r_net, k_diss, C_sat, C_bulk=0.0):
    """Cobertura superficial en estado estacionario (metal limpio por defecto)."""
    denom = k_diss * (1.0 - C_bulk / C_sat)
    if denom <= 0:
        return math.inf
    return r_net / denom


def t_saturation_h(r_net, C_sat, rho=None):
    """Tiempo caracteristico de saturacion [h] = C_sat/(r_net*rho)."""
    if rho is None:
        rho = RHO_ACUM
    if r_net <= 0:
        return math.inf
    return (C_sat / (r_net * rho)) / 3600.0


def integrate_accumulation(r_net, C_sat, rho=None, dt_h=1.0,
                           t_max_h=None, sat_frac=0.999):
    """Integra dC/dt = r_net*rho*(1-C/C_sat) (Euler, paso dt_h h).
    Acota (1-C/C_sat) a [0,1]; si C_bulk supera C_sat se detiene.
    Devuelve (t_h, C_bulk[mol/m3], CF)."""
    if rho is None:
        rho = RHO_ACUM
    if r_net <= 0:
        out = ([0.0], [0.0], [1.0])
        return tuple(np.array(x) for x in out) if _HAS_NP else out

    tsat = t_saturation_h(r_net, C_sat, rho)
    if tsat < 5 * dt_h:
        print("[AVISO] t_sat = %.4f h < 5*dt (%.1f h). Perfil Euler impreciso; "
              "reduzca dt_h (p. ej. dt_h=%.5f)." % (tsat, 5 * dt_h, tsat / 20.0))
    if t_max_h is None:
        t_max_h = 8.0 * tsat

    dt_s = dt_h * 3600.0
    incr = r_net * rho * dt_s
    t_list, C_list, CF_list = [0.0], [0.0], [1.0]
    C, t, steps = 0.0, 0.0, 0
    while C < sat_frac * C_sat and t < t_max_h and steps < 5_000_000:
        cf = max(0.0, min(1.0, 1.0 - C / C_sat))
        if cf <= 0.0:
            break
        C = min(C + incr * cf, C_sat)
        t += dt_h
        steps += 1
        t_list.append(t); C_list.append(C)
        CF_list.append(max(0.0, min(1.0, 1.0 - C / C_sat)))
    if _HAS_NP:
        return np.array(t_list), np.array(C_list), np.array(CF_list)
    return t_list, C_list, CF_list


# =============================================================================
# 6. ENTRADAS: conversion y lectura de CSV
# =============================================================================
def p_from_mole_fractions(y, P_total=P_BASE):
    return {sp: y[sp] * P_total for sp in y}


def p_from_molar_flows(F, P_total=P_BASE):
    Ftot = sum(F.values())
    return {sp: (F[sp] / Ftot) * P_total for sp in F}


def read_profile_csv(path, colmap, kind="pressure", P_total=P_BASE, delimiter=","):
    """Lee un perfil axial de un CSV de DWSIM.
    colmap: {especie: nombre_columna}. kind: 'pressure'|'fraction'|'flow'."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter=delimiter):
            raw = {sp: float(row[col]) for sp, col in colmap.items()}
            if kind == "pressure":
                rows.append(raw)
            elif kind == "fraction":
                rows.append(p_from_mole_fractions(raw, P_total))
            elif kind == "flow":
                rows.append(p_from_molar_flows(raw, P_total))
            else:
                raise ValueError("kind: 'pressure', 'fraction' o 'flow'")
    return rows


# =============================================================================
# 7. CASO UNICO Y BARRIDO EN TEMPERATURA
# =============================================================================
def run_case(profile, T=T_BASE, scenarios=None, dt_h=1.0, mode="rates", verbose=True):
    """Caso a una sola temperatura. profile: dict (un punto) o lista (perfil)."""
    if scenarios is None:
        scenarios = SCENARIOS
    K = kinetic_constants(T)
    rates = mean_carbon_rates(profile, T, K, mode=mode)
    out = {"T": T, "rates": rates, "scenarios": {}}
    if verbose:
        print("=" * 66)
        print("  T = %.2f K   |  r_crack=%+.3e  r_Boud=%+.3e  r_net=%+.3e"
              % (T, rates["r_crack"], rates["r_Boud"], rates["r_net"]))
        if rates["r_Boud"] <= 0:
            print("  -> Boudouard no aporta carbono a esta T.")
    for name, sc in scenarios.items():
        th = theta_ss(rates["r_net"], sc["k_diss"], sc["C_sat"])
        tsat = t_saturation_h(rates["r_net"], sc["C_sat"])
        t_h, C_b, CF = integrate_accumulation(rates["r_net"], sc["C_sat"], dt_h=dt_h)
        out["scenarios"][name] = dict(k_diss=sc["k_diss"], C_sat=sc["C_sat"],
            theta_ss=th, t_sat_h=tsat, t_h=t_h, C_bulk=C_b, CF=CF)
        if verbose:
            warn = "  [theta_SS>1: no fisico]" if th > 1 else ""
            print("  [%-11s] theta_SS=%.3e  t_sat=%.4g h%s"
                  % (name, th, tsat, warn))
    if verbose:
        print("=" * 66)
    return out


def run_sweep(data, scenarios=None, dt_h=1.0, mode="rates"):
    """Barrido en T. data: dict {T_K: profile} (profile = dict o lista de p [Pa])."""
    return {T: run_case(data[T], T=T, scenarios=scenarios, dt_h=dt_h,
                        mode=mode, verbose=False) for T in sorted(data)}


# =============================================================================
# 8. SALIDA: tabla, CSV y graficas
# =============================================================================
def print_sweep(sweep, scenario="central"):
    print("=" * 84)
    print(" Barrido en T   (escenario '%s')" % scenario)
    print(" %6s %13s %13s %13s %11s %12s"
          % ("T[K]", "r_crack", "r_Boud", "r_net", "theta_SS", "t_sat[h]"))
    for T in sorted(sweep):
        r = sweep[T]["rates"]; s = sweep[T]["scenarios"][scenario]
        print(" %6.0f %13.3e %13.3e %13.3e %11.2e %12.4g"
              % (T, r["r_crack"], r["r_Boud"], r["r_net"], s["theta_ss"], s["t_sat_h"]))
    print("=" * 84)


def export_sweep_csv(path, sweep, scenarios=("conservador", "central", "optimista")):
    """CSV resumen: una fila por temperatura, columnas por escenario."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["T (K)", "r_crack (mol/kg/s)", "r_Boud (mol/kg/s)", "r_net (mol/kg/s)"]
        for sc in scenarios:
            header += ["theta_SS_%s" % sc, "t_sat_%s (h)" % sc]
        w.writerow(header)
        for T in sorted(sweep):
            r = sweep[T]["rates"]
            row = [T, r["r_crack"], r["r_Boud"], r["r_net"]]
            for sc in scenarios:
                s = sweep[T]["scenarios"][sc]
                row += [s["theta_ss"], s["t_sat_h"]]
            w.writerow(row)


def export_accumulation_csv(path, result, scenario="central"):
    """CSV del perfil temporal [t(h), C_bulk, CF] para un caso/escenario."""
    s = result["scenarios"][scenario]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t (h)", "C_bulk (mol/m3)", "CF"])
        for row in zip(s["t_h"], s["C_bulk"], s["CF"]):
            w.writerow(list(row))


def plot_sweep(sweep, scenario="central", path=None):
    """t_sat vs T (requiere matplotlib)."""
    if not _HAS_MPL:
        print("[AVISO] matplotlib no disponible; se omite la grafica.")
        return
    Ts = sorted(sweep)
    tsat = [sweep[T]["scenarios"][scenario]["t_sat_h"] for T in Ts]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(Ts, tsat, "o-")
    ax.set_xlabel("T (K)"); ax.set_ylabel("t_sat (h)")
    ax.set_title("Tiempo de saturacion vs T (%s)" % scenario)
    fig.tight_layout()
    fig.savefig(path, dpi=150) if path else plt.show()


# =============================================================================
# 8b. LECTOR DEL FORMATO DWSIM (.xlsx con perfil axial por temperatura)
# =============================================================================
# Mapeo de especies a los nombres de columna que exporta DWSIM.
DWSIM_COLMAP = {
    "CH4": "Methane MolFrac",
    "CO2": "Carbon dioxide MolFrac",
    "CO":  "Carbon monoxide MolFrac",
    "H2":  "Hydrogen MolFrac",
    "H2O": "Water MolFrac",
}


def read_dwsim_sheet(path, sheet=0, L_max=H_LM, colmap=None):
    """Lee una hoja de DWSIM (perfil axial a una temperatura).

    Usa las fracciones molares y la presion total: p_i = y_i * P.
    Trunca a la zona reactiva (Length <= L_max = altura de Sn, 0.6 m).
    Devuelve (T [K], profile) donde profile es lista de dicts de p [Pa].
    """
    if not _HAS_PD:
        raise ImportError("Se necesita pandas para leer .xlsx: "
                          "pip install pandas openpyxl")
    if colmap is None:
        colmap = DWSIM_COLMAP
    df = pd.read_excel(path, sheet_name=sheet)
    T = float(df["Temperature (K)"].iloc[len(df) // 2])  # T constante de la hoja
    df = df[df["Length (m)"] <= L_max + 1e-9]
    profile = []
    for _, row in df.iterrows():
        P = float(row["Pressure (Pa)"])
        profile.append({sp: float(row[col]) * P for sp, col in colmap.items()})
    return T, profile


def load_dwsim_files(paths, L_max=H_LM):
    """Una corrida (= una T) por fichero .xlsx. paths: lista de rutas.
    Devuelve data = {T: profile} listo para run_sweep()."""
    data = {}
    for p in paths:
        T, prof = read_dwsim_sheet(p, L_max=L_max)
        data[round(T)] = prof
    return data


def load_dwsim_workbook(path, L_max=H_LM):
    """Una corrida (= una T) por HOJA dentro de un mismo .xlsx.
    Devuelve data = {T: profile} listo para run_sweep()."""
    if not _HAS_PD:
        raise ImportError("Se necesita pandas para leer .xlsx.")
    xls = pd.ExcelFile(path)
    data = {}
    for sh in xls.sheet_names:
        T, prof = read_dwsim_sheet(path, sheet=sh, L_max=L_max)
        data[round(T)] = prof
    return data


# =============================================================================
# 9. PUNTO DE ENTRADA  ---  CARGA LOS DATOS DE DWSIM
# =============================================================================
if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # OPCION A: una corrida de DWSIM por fichero .xlsx (una por temperatura)
    #   Pon aqui las rutas de tus 12 ficheros (723..1273 K):
    DWSIM_FILES = [
        # "Datos_DWSIM_723.xlsx",
        # "Datos_DWSIM_773.xlsx",
        # ...
        # "Datos_DWSIM_1273.xlsx",
    ]
    # OPCION B: todas las temperaturas como hojas de un mismo .xlsx:
    #   data = load_dwsim_workbook("Datos_DWSIM_barrido.xlsx")
    # ------------------------------------------------------------------ #

    if DWSIM_FILES:
        data = load_dwsim_files(DWSIM_FILES)        # {T: profile}
    else:
        # Demo con el fichero de ejemplo (una sola T = 873 K):
        try:
            T, prof = read_dwsim_sheet("Ejemplo_Datos_DWSIM.xlsx")
            data = {round(T): prof}
        except Exception:
            data = {}

    if not data:
        print("Indica tus ficheros de DWSIM en DWSIM_FILES (o usa "
              "load_dwsim_workbook).")
    else:
        # mode='rates' (riguroso) o 'pressures' (equivalente al Excel)
        sweep = run_sweep(data, dt_h=1.0, mode="rates")
        print_sweep(sweep, scenario="central")
        export_sweep_csv("carbono_barrido.csv", sweep)
        print("Resumen del barrido -> carbono_barrido.csv")
        Tref = 1173 if 1173 in sweep else max(sweep)
        export_accumulation_csv("carbono_acumulacion_%dK.csv" % Tref, sweep[Tref])
        print("Acumulacion C_bulk(t) a %d K -> carbono_acumulacion_%dK.csv"
              % (Tref, Tref))
        # plot_sweep(sweep, "central", "tsat_vs_T.png")
