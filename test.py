Vao = 1
Va = 1
gamma = 1.4
k0 = 0.7
mb = 0.15           # Masse fusée (kg)
A = 0.00785        # Surface de la fusée (m²)
Ae = 0.00005024    # Surface de la buse (m²)
rho_w = 1000       # Densité de l'eau (kg/m³)
patm = 101325      # Pression atmosphérique (Pa)
rho_atm = 1.225    # Densité de l'air (kg/m³)

# 1. CONSTANTES PHYSIQUES
g = 9.81               # Gravité (m/s²)
rho_w = 1.0e3         # Densité eau (kg/m³)
Cd = 0.35             # Coefficient traînée

# 2. PARAMÈTRES GÉOMÉTRIQUES
D = 0.1               # Diamètre fusée (m)
De = 0.008            # Diamètre buse (m)
V = 0.0015           # Volume total (m³)
p_ino = 500000        # Pression initiale (Pa)

# 3. FONCTIONS DE CALCUL DES FORCES
def Weight(mw):
    F_weight = (mb + mw) * g

def Drag(v):
    FD = 0.5 * rho_atm * v**2 * Cd * A

def Thrust(v_e):
    F_thrust = rho_w * Ae * v_e**2

# 4. PRESSION INTERNE
def internal_pressure(Vw):
    p_in = p_ino * (Vao / Va)**gamma

# 5. VITESSE DE SORTIE DE L'EAU
def exit_velocity(Vw):
    k = Vw / V
    v_e = ((2*(p_ino*((1-k0)/(1-k))**(gamma)-patm))/((rho_w)*(1-((Ae)/(A))**2)))**(1/2) 