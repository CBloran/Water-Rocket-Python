
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import json
# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.81               # Accélération de la gravité (m/s²)
rho_w = 1.0e3         # Densité de l'eau (kg/m³)
rho_atm = 1.23        # Densité de l'air (kg/m³)
Cd = 0.831             # Coefficient de traînée aérodynamique
gamma = 1.4           # Coefficient adiabatique de l'air
patm = 101325         # Pression atmosphérique (Pa)
R_air = 287.05

# Paramètres géométriques et initiaux du rocket
#les noms avec o sont des paramètres à t = 0s (on peut pas mettre 0 dans la variable donc o = initial, tandis que "in" comme dans p_ino est pour "intérieur" ou "interne")
# les paramètres qu'on ne connait pas sont notés float pour l'instant
D = float(0.09)               # Diamètre de la fusée (m)
De = float(0.008)            # Diamètre de la buse (m)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
V = float(0.0015)          # Volume total du rocket (m³)
#p_ino = float(500000)           # Pression initiale à l'intérieur (Pa)
T_initial = 293.15  # 20°C
#p_ino += 100000             # Convertir la pression relative en pression absolue
mb = float(0.3)             # Masse structurelle (kg)
#Vw0 = float(0.0009)              # Volume d'eau initial dans la fusée

#=============================================================
#--------------------MODEL CORECTION--------------------------
#=============================================================
ForceOffset = 0#-7

# ============================================================
#  FORCES EN PRESENCE
# ============================================================
def Weight(mw):
    """
    mw : water mass
    Calculate the weight of the rocket
    """
    
    W = (mb + mw) * g
    return W
def Drag(v):
    """
    v : Rocket speed
    Calcculate the drag force
    """
    FD = 0.5 * rho_atm * v**2 * Cd * A
    return FD
def Thrust(v_e):
    """
    v_e : exit velocity of water
    calculate the thrust force based on the exit velocity of water
    """
    F_thrust = rho_w * Ae * v_e**2 
    F_thrust = max(F_thrust+ForceOffset, 0)
    return F_thrust
#===============================================================
def water_exit_velocity(_p_in):
    """
    k : ratio of remaining water volume to total volume
    p_in : internal pressure
    """
    try:
        v_e = ((2*(_p_in-patm))/((rho_w)*(1-(Ae/A)**2)))**(1/2) #calculate the exit velocity of water based on bernouilli's equation
        #  OTHER BERNOUILLI EQUATION  #
        #delta_P = p_in - patm
        #v_e = (2.0 * delta_P / (rho_w*(1-(Ae/A)**2)))**(1/2) #calculate the exit velocity of water based on bernouilli's equation
        print(v_e)
        return v_e
    except:
        return 0

def equation_water_mass(v_e):
    """
    v_e : water ejection velocity
    Calculate the mass of water in the rocket based on the internal pressure
    """
    dmw_dt = -rho_w * Ae * v_e
    return dmw_dt

def equation_vel_air(v, p_in, _p_ino):
    F_drag = Drag(v)
    F_weight = Weight(0)
    
    if p_in > patm:
        # Température de l'air restant
        T_air = T_initial * (p_in / _p_ino) ** ((gamma - 1) / gamma)
        
        # Rapport de pression critique pour le col sonique
        pressure_ratio_critical = ((gamma + 1) / 2) ** (gamma / (gamma - 1))
        
        if p_in / patm >= pressure_ratio_critical:
            # Écoulement sonique (bloqué)
            v_e = math.sqrt(gamma * R_air * T_air)
            rho_e = (p_in / (R_air * T_air)) * (2 / (gamma + 1)) ** (1 / (gamma - 1))
            F_thrust = 0#Ae * (rho_e * v_e**2 + (p_in * (2/(gamma+1))**(gamma/(gamma-1)) - patm))
        else:
            # Écoulement subsonique
            v_e = math.sqrt(2 * gamma * R_air * T_air / (gamma - 1) * 
                          (1 - (patm / p_in) ** ((gamma - 1) / gamma)))
            F_thrust = 0#Ae * v_e**2  # Approximation
    else:
        F_thrust = 0
        v_e = 0
    
    
    return F_thrust, v_e
    
#=======================================================
########           Runguert-Kutta             ##########
#=======================================================

def systeme_complet(t, y, _p_ino, deltaT, mw0):
        """
        Params:
        t = time
        y = [h, v, mw, p_in, F]
        return [dh_dt, dv_dt, dv_e_dt, dmw_dt, dp_in_dt, F]
        """
        
        def internal_pressure(mw):
            """
            Vw : volume d'eau restant
            Calcul de la pression interne en fonction du volume d'eau restant.
            """
            Va = V - mw/rho_w  # remaining air volume
            Vao = V - mw0/rho_w  # initial air volume

            if Va > 0:
                p_in = _p_ino * (Vao / Va)**gamma # calculate the remaining internal pressure based on adiabatic law
                return p_in
            else:
                return 0 
            
        h, v, mw, v_e, p_in, F = y

        NextP_in = internal_pressure(mw)
        dm_dt = equation_water_mass(v_e)

        F_drag = Drag(v)
        F_weight = Weight(mw)

              
        if mw > 0.00000001:
                Nextv_e = water_exit_velocity(p_in)
                
                F_thrust = Thrust(v_e)
                F = F_thrust - F_weight
                dv_dt = (F-F_drag)/(mb + mw)

        elif p_in > patm:
            
            F_thrust, Nextv_e = equation_vel_air(v, p_in, _p_ino)
            
            dv_dt = (F_thrust-F_weight -F_drag)/(mb + mw)
            F = F_thrust - F_weight
        else:
            Nextv_e = 0
            dv_dt = (-F_weight -F_drag)/(mb + mw)
            F = -F_weight

        dh_dt = v    




        return [dh_dt, dv_dt, dm_dt, Nextv_e, NextP_in, F]




def Rungekutta(_mw0, _p_ino, stepNbr=50000):
    """
    stepNbr : number of points in the simulation
    """

    # Initial conditions
    t0 = 0.00001 # Initial time
    h0 = 0.00001 # Initial height
    v0 = 0 # Initial velocity
    mw0 = _mw0 # Initial water volume
    p_ino = _p_ino + patm # Initial pressure

    
    
    
    
    # Initialisation
    ts = [t0]
    Fs = [0]
    ys = [[h0, v0, mw0, 0, p_ino, 0]]  # Stocker toutes les variables dans une liste
    
    deltaT = 0.001  # Pas de temps constant
    i = 0
    
    while ys[-1][0] > 0 and i < stepNbr :
        i += 1
        t_current = ts[-1]
        y_current = ys[-1]
        #print(y_current[0])
        
        
        #print(p_in)
        # RK4 standard
        k1 = systeme_complet(t_current, y_current, p_ino, deltaT, mw0)
        
        k2 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k1[j] for j in range(4)]+[y_current[4]]+[k1[4]]+[k1[5]], p_ino, deltaT/2, Vw0)
        
        k3 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k2[j] for j in range(4)]+[y_current[4]]+[k1[4]]+[k1[5]], p_ino, deltaT/2, Vw0)
        
        k4 = systeme_complet(t_current + deltaT, [y_current[j] + deltaT * k3[j] for j in range(4)]+[y_current[4]]+[k1[4]]+[k1[5]], p_ino, deltaT, Vw0)
        
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(4)]+[k1[4]] + [k1[5]]
            
        ts.append(t_current + deltaT)
        ys.append(y_new)


    
    
    #graphMathPlot(ts, ys)
    # Spliting the results
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    v_es = [y[2] for y in ys]
    mws = [y[3] for y in ys]
    ps = [y[4] for y in ys]
    Fs = [y[5] for y in ys]

    with open("height.json", "w", encoding="utf-8") as f :
        for y in ys :
            strin = str(y)
            f.write(strin + "\n")
    with open("time", "w", encoding="utf-8") as f:
        for t in ts :
            strin = str(t)
            f.write(strin + "\n")

    #print(Fs)
    return ts, ys


def RungekuttaMax(x, y):
    result = Rungekutta(x, y)[1]
    hs = [y[0] for y in result]
    return max(hs)



#print(tsys[2])
#txtGraph(tsys[0], tsys[1], "Height.txt")

#=======================================================
########           Graph V2 MathPlot             #########
#=======================================================
def graphMathPlot(Times, dataIn, additionalValues=None):
    

    data = dataIn

    x = Times 

    plt.figure(figsize=(10, 5))

    # Séparation des composantes
    h  = [row[0] for row in data]
    plt.plot(x, h,  label="h",  color="red",   marker="o", markersize = 2)

    if len(data[0]) > 1:
        v  = [row[1] for row in data]
        plt.plot(x, v,  label="v",  color="blue",  marker="o", markersize = 2)

    if len(data[0]) > 2:
        v_e = [row[2]*1000 for row in data]
        plt.plot(x, v_e, label="v_e", color="green", marker="o", markersize = 2)

    if len(data[0]) > 3:
        mw = [row[3]/1000 for row in data]
        plt.plot(x, mw, label="mw", color="purple", marker="o", markersize = 2)
    
    if len(data[0]) > 4:
        ps = [row[4] for row in data]
        plt.plot(x, ps, label="p_in", color="orange", marker="o", markersize = 2)
    if len(data[0]) > 5:
        F = [row[5] for row in data]
        plt.plot(x, F, label="F", color="purple", marker="o", markersize = 2)
    

    
    
    
    
    if additionalValues:
        plt.plot(x, additionalValues, label="F_thrust", color="orange", marker="o", markersize = 2)

    plt.xlabel("Index")
    plt.ylabel("Valeur")
    plt.title("Évolution de h, v et Vw")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()
#====================================================================================
#-----------------------------------GRAPHIC--3D--------------------------------------
#====================================================================================
def graphique_3d_parametres(fonction, x_min, x_max, y_min, y_max, nb_points_x=50, nb_points_y=50):
    """
    Crée un graphique 3D d'une fonction à deux paramètres.
    
    Paramètres:
    -----------
    fonction : function
        La fonction à évaluer. Doit prendre deux arguments (x, y) et retourner une valeur.
    x_min, x_max : float
        Bornes de l'intervalle pour le paramètre x
    y_min, y_max : float
        Bornes de l'intervalle pour le paramètre y
    nb_points_x, nb_points_y : int
        Nombre de points à générer sur chaque axe (défaut: 50)
    """
    
    # Création des grilles de points
    maxVolume = 0
    maxHeight = 0
    
    x = np.linspace(x_min, x_max, nb_points_x)
    y = np.linspace(y_min, y_max, nb_points_y)
    
    # Création de la grille 2D
    X, Y = np.meshgrid(x, y)
    
    # Calcul des valeurs Z point par point
    Z = np.zeros_like(X)
    for i in range(nb_points_y):
        for j in range(nb_points_x):
            Z[j, i] = fonction(X[j, i], Y[j, i])
            if Z[j, i] > maxHeight:
                maxHeight = Z[j, i]
                maxVolume = X[j, i]
            else:
                print(maxVolume)
        
        
    # Création de la figure 3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Tracé de la surface
    surface = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    
    # Personnalisation du graphique
    ax.set_xlabel("Volume d'eau (m³)")
    ax.set_ylabel("Pression d'entrée (Pa)")
    ax.set_zlabel('Hauteur atteinte')
    ax.set_title("Evolution de la hauteur atteinte en fonction du volume d'eau et de la pression")
    
    # Ajout d'une barre de couleur
    plt.colorbar(surface, ax=ax, shrink=0.5, aspect=5)
    
    plt.show()
    
    return fig, ax
#graphique_3d_parametres(RungekuttaMax, 0.0001, 0.0010, 100000, 500000)

def batch():
    array = []
    i = 100000
    while i <= 500000:
        print(i)
        hmax = RungekuttaMax(0.0005, i)
        array.append([i, hmax])
        i += 1000
    with open("AllPressure.txt", "w", encoding="utf-8") as f:
        for y in array:
            strin = str(y)
            strin = strin.replace("[", "")
            strin = strin.replace("]", "")
            
            f.write(strin + "\n")

def FinBestPressure(Height):
    with open("AllPressure.txt", "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                array = line.split(",")
                if float(array[1]) >= Height:
                    return float(array[0])

if __name__ == "__main__":
    print("Veuillez choisir le mode d'execution :")
    print("1. Mode Volume d'eau constant")
    print("2. Mode graphique 3D")
    print("3. Géneration du résumer pression - hauteur")
    print("4. Trouver la meilleur pression pour la hauteur")
    choice = int(input())
    if choice == 1:
        print("Veuillez entrer le masse d'eau (en kg) :")
        Vw0 = float(input())
        print("Veuillez entrer la pression d'entrée (en Pa) :")
        p_ino = float(input())
        tsys = Rungekutta(Vw0, p_ino)
        print("Hauteur maximale atteinte : ", max(tsys[1]))
        graphMathPlot(tsys[0], tsys[1])
    elif choice == 2:
        graphique_3d_parametres(RungekuttaMax, 0.0001, 0.0010, 100000, 500000)
    elif choice == 3:
        batch()
    elif choice == 4:
        print("Veuillez entrer la hauteur maximale souhaitée (en m) :")
        Height = float(input())
        print("La pression d'entrée optimale pour atteindre une hauteur de ", Height, "m est de : ", FinBestPressure(Height), "Pa")