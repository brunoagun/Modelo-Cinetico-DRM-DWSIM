import sys; sys.path.insert(0,'.')
import pandas as pd, math
import carbono_caso_base as m

PATH = '/mnt/user-data/uploads/ResultadosDWSIM_SimPlevan.xlsx'
H_LM = m.H_LM  # 0.6 m

def find_col(cols, key, kind):
    for c in cols:
        cl = c.lower().replace(' ','')
        if key in cl and kind in cl:
            return c
    return None

SPEC = {'CH4':'methane','CO2':'dioxide','CO':'monoxide','H2':'hydrogen','H2O':'water'}

def read_sheet(path, sheet):
    df = pd.read_excel(path, sheet_name=sheet)
    cols = list(df.columns)
    frac = {sp: find_col(cols, k, 'molfrac') for sp,k in SPEC.items()}
    flow = {sp: find_col(cols, k, 'molflow') for sp,k in SPEC.items()}
    T = float(df['Temperature (K)'].iloc[len(df)//2])
    Lcol = 'Length (m)'; Pcol = 'Pressure (Pa)'
    df = df[df[Lcol] <= H_LM + 1e-9].reset_index(drop=True)
    return df, T, frac, flow, Lcol, Pcol

xls = pd.ExcelFile(PATH)
sheets = [s for s in xls.sheet_names if s.upper().startswith('PFR')]

rows = []
axial = {}
for sh in sheets:
    df, T, frac, flow, Lcol, Pcol = read_sheet(PATH, sh)
    K = m.kinetic_constants(T)
    # perfil de presiones parciales y velocidades
    prof, rc, rb, rn, aff_c, aff_b, aC, X_CH4, X_CO2 = [],[],[],[],[],[],[],[],[]
    L = df[Lcol].values
    F_CH4_0 = df[flow['CH4']].iloc[0]; F_CO2_0 = df[flow['CO2']].iloc[0]
    for i in range(len(df)):
        P = float(df[Pcol].iloc[i])
        p = {sp: float(df[frac[sp]].iloc[i])*P for sp in SPEC}
        prof.append(p)
        r = m.carbon_rates(p, T, K)
        rc.append(r['r_crack']); rb.append(r['r_Boud']); rn.append(r['r_net'])
        aff_c.append(r['aff_crack']); aff_b.append(r['aff_boud'])
        # actividad de carbono via craqueo: a_C = p_H2^2/(Keq_crack*p_CH4) [bar]
        bH2=p['H2']/1e5; bCH4=p['CH4']/1e5
        aC.append((bH2**2)/(K['Keq_crack']*bCH4) if bCH4>0 else float('nan'))
        X_CH4.append(1 - df[flow['CH4']].iloc[i]/F_CH4_0)
        X_CO2.append(1 - df[flow['CO2']].iloc[i]/F_CO2_0)
    # promedios de lecho
    rates_avg = m.mean_carbon_rates(prof, T, K, mode='rates')
    # caracterizacion axial
    rn_arr = rn
    i_peak = max(range(len(rn_arr)), key=lambda i: rn_arr[i])
    # longitud para alcanzar 99% de la conversion final de CH4
    Xf = X_CH4[-1]
    L99 = next((L[i] for i in range(len(L)) if Xf>0 and X_CH4[i]>=0.99*Xf), L[-1])
    # outlet
    out = prof[-1]
    H2CO = out['H2']/out['CO'] if out['CO']>0 else float('inf')
    # escenarios
    sc = {}
    for name,s in m.SCENARIOS.items():
        sc[name] = dict(theta=m.theta_ss(rates_avg['r_net'], s['k_diss'], s['C_sat']),
                        tsat=m.t_saturation_h(rates_avg['r_net'], s['C_sat']))
    rows.append(dict(T=round(T,2), Xch4=X_CH4[-1], Xco2=X_CO2[-1], H2CO=H2CO,
        p_CH4=out['CH4'], p_CO2=out['CO2'], p_CO=out['CO'], p_H2=out['H2'], p_H2O=out['H2O'],
        r_inlet=rn[0], r_outlet=rn[-1], r_avg=rates_avg['r_net'],
        rc_avg=rates_avg['r_crack'], rb_avg=rates_avg['r_Boud'],
        aff_in=aff_c[0], aff_out=aff_c[-1], aC_out=aC[-1], affb_out=aff_b[-1],
        L_peak=L[i_peak], L99=L99,
        th_cons=sc['conservador']['theta'], th_cent=sc['central']['theta'], th_opt=sc['optimista']['theta'],
        ts_cons=sc['conservador']['tsat'], ts_cent=sc['central']['tsat'], ts_opt=sc['optimista']['tsat']))
    axial[round(T)] = dict(L=list(L), X=X_CH4, rn=rn, rc=rc, rb=rb, aff=aff_c, aC=aC)

res = pd.DataFrame(rows).sort_values('T').reset_index(drop=True)
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 50)

print('='*120)
print('RESUMEN POR TEMPERATURA (promedios de lecho, zona reactiva L<=0.6 m)')
print('='*120)
show = res[['T','Xch4','Xco2','H2CO','r_inlet','r_outlet','r_avg','rb_avg','aff_in','aff_out','aC_out','ts_cent']].copy()
show.columns = ['T[K]','X_CH4','X_CO2','H2/CO','r_in','r_out','r_avg','rB_avg','aff_in','aff_out','aC_out','tsat_c[h]']
print(show.to_string(index=False, float_format=lambda x: '%.4g'%x))
print()
print('Presiones parciales de SALIDA [Pa]:')
print(res[['T','p_CH4','p_CO2','p_CO','p_H2','p_H2O']].to_string(index=False, float_format=lambda x:'%.0f'%x))
print()
print('t_sat [h] por escenario:')
print(res[['T','ts_cons','ts_cent','ts_opt']].to_string(index=False, float_format=lambda x:'%.4g'%x))
print()
print('theta_SS por escenario (chequeo <=1):')
print(res[['T','th_cons','th_cent','th_opt']].to_string(index=False, float_format=lambda x:'%.3e'%x))
print()
print('Localizacion del pico de r_net y L99 (longitud al 99%% de conversion):')
print(res[['T','L_peak','L99']].to_string(index=False, float_format=lambda x:'%.3f'%x))

res.to_csv('analisis_resumen.csv', index=False)
# guardar perfiles axiales
import json
with open('analisis_axial.json','w') as f: json.dump(axial, f)
print('\nGuardado: analisis_resumen.csv, analisis_axial.json')
