import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques
g = 9.81               # m/s²
rho_eau = 1000.0       # kg/m³
rho_air = 1.225        # kg/m³ (air à 15°C, 1 atm)
p_atm = 101325.0       # Pa
gamma = 1.4            # rapport des chaleurs spécifiques de l'air
n_poly = 1.3           # exposant polytropique (entre 1.2 et 1.4)

# Paramètres par défaut de la fusée (à modifier selon votre modèle)
M_vide = 0.100         # kg (masse de la structure vide, sans eau ni air)
V_total = 1.5e-3       # m³ (1.5 L)
diametre = 0.09        # m (diamètre du corps)
diametre_tuyere = 0.009 # m (9 mm)
C_d = 0.6              # coefficient de traînée (estimé)

# Paramètres initiaux (modifiables)
p_init_bar = 5.0       # pression initiale en bar
ratio_eau = 0.5        # fraction du volume occupée par l'eau (ex: 30%)

# Conversion
p_init = p_init_bar * 1e5  # Pa
V_eau0 = V_total * ratio_eau
V_air0 = V_total - V_eau0
masse_air0 = (p_init * V_air0) / (287.0 * 293.0)  # kg (loi des gaz parfaits, T=20°C)

# Surfaces
A_ref = np.pi * (diametre/2)**2          # section de référence pour la traînée
A_tuyere = np.pi * (diametre_tuyere/2)**2

# Conditions initiales pour la phase propulsée
t = 0.0
h = 0.0
v = 0.0
V_eau = V_eau0
p = p_init
masse_air = masse_air0   # on suppose que l'air ne s'échappe pas pendant la phase eau
masse_tot = M_vide + masse_air + rho_eau * V_eau

# Listes pour enregistrer les données (optionnel)
temps = []
altitude = []
vitesse = []

# Pas de temps pour Euler
dt = 0.0005  # 0.5 ms (ajuster si nécessaire)

print("Début de la phase de propulsion par eau...")

while V_eau > 0 and p > p_atm:
    # Mise à jour de la masse et de la pression
    V_air = V_total - V_eau
    p = p_init * (V_air0 / V_air)**n_poly   # loi polytropique
    if p <= p_atm:
        p = p_atm
        v_e = 0.0
    else:
        # Vitesse d'éjection (Bernoulli)
        v_e = np.sqrt(2 * (p - p_atm) / rho_eau)
    
    # Dérivées
    dV_eau = - A_tuyere * v_e
    F = rho_eau * A_tuyere * v_e**2
    D = 0.5 * rho_air * C_d * A_ref * v**2
    dv = (F - D - masse_tot * g) / masse_tot
    dh = v
    
    # Euler
    V_eau += dV_eau * dt
    if V_eau < 0:
        V_eau = 0.0
    v += dv * dt
    h += dh * dt
    t += dt
    
    # Mise à jour de la masse
    masse_tot = M_vide + masse_air + rho_eau * V_eau
    
    # Enregistrement
    temps.append(t)
    altitude.append(h)
    vitesse.append(v)

# Fin de la phase propulsée
v_b = v
h_b = h
t_b = t
print(f"Fin de la phase propulsée à t = {t_b:.3f} s")
print(f"Vitesse en fin de propulsion : {v_b:.2f} m/s")
print(f"Altitude en fin de propulsion : {h_b:.2f} m")

# Phase balistique (avec traînée) – formule analytique
M_final = M_vide + masse_air   # on garde la masse d'air (négligeable)
psi = (rho_air * C_d * A_ref) / (2 * M_final)   # ψ = D/M
if v_b > 0:
    # Formule (35) du document Fischer et al.
    h_bal = (1/(2*psi)) * np.log(1 + (psi/g) * v_b**2)
else:
    h_bal = 0.0

h_max = h_b + h_bal
print(f"Altitude gagnée en phase balistique : {h_bal:.2f} m")
print(f"Altitude maximale atteinte : {h_max:.2f} m")

# Tracé optionnel
plt.figure(figsize=(8,5))
plt.plot(temps, altitude, label='Altitude (m)')
plt.xlabel('Temps (s)')
plt.ylabel('Altitude (m)')
plt.title('Simulation de vol d\'une fusée à eau')
plt.grid(True)
plt.legend()
plt.show()