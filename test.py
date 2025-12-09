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
F_thrust = 0
F_drag = 0
F_weight = 0
mw_current = 0
Vw = 0.0005          # Volume d'eau initial (m³)

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


def systeme_complet(v, v_e):
    
    dh_dt = v
    dVw_dt = -Ae * ((2*(p_ino*((1-k0)/(1-k))**(gamma)-patm))/((rho_w)*(1-((Ae)/(A))**2)))**(1/2) 
    dv_dt = (F_thrust - F_drag - F_weight) / (mb + mw_current)


    return [dh_dt, dv_dt, dVw_dt]

t_current = ts[-1]
y_current = ys[-1]
deltaT = 0.01  # Pas de temps constant
i = 0
def ronguert_kutta(stepNbr=3000):
    # Initial conditions
    t0 = 0.00001 # Initial time
    h0 = 0.00001 # Initial height
    v0 = 0 # Initial velocity
    Vw0 = 0.0005
    
    # Initialisation
    ts = [t0]
    ys = [[h0, v0, Vw0]]  # Stocker toutes les variables dans une liste
    
    while ys[-1][0] > 0 and i < stepNbr:
        t_current = ts[-1]
        y_current = ys[-1]
        
        # RK4 standard
        k1 = systeme_complet(t_current, y_current)
        k2 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k1[j] for j in range(3)])
        k3 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k2[j] for j in range(3)])
        k4 = systeme_complet(t_current + deltaT, [y_current[j] + deltaT * k3[j] for j in range(3)])
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(3)]
        ts.append(t_current + deltaT)
        ys.append(y_new)
    ...











    