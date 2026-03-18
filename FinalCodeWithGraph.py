
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
# ============================================================
def internal_pressure(_p_ino, Vw, Vw0):
    """
    Vw : volume d'eau restant
    Calcul de la pression interne en fonction du volume d'eau restant.
    """
    Va = V - Vw  # remaining air volume
    Vao = V - Vw0  # initial air volume

    if Va > 0:
        p_in = _p_ino * (Vao / Va)**gamma # calculate the remaining internal pressure based on adiabatic law
        return p_in
    else:
        return 0 

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

def equation_vel(v, Vw, Vw0, _p_ino):
            
            mw_current = rho_w * Vw  # Masse d'eau actuelle
            F_drag = Drag(v)
            F_weight = Weight(mw_current)  
            p_in = internal_pressure(_p_ino, Vw, Vw0)
            k0 = Vw0 / V

            if Vw > 0:
                
                v_e = water_exit_velocity(p_in)
                F_thrust = Thrust(v_e)
                #print(F_thrust)
            else:
                F_thrust = 0
            
            if Vw > 0:
                dv_dt = (F_thrust - F_drag - F_weight) / (mb + mw_current)
            else:
                dv_dt = (-F_drag - F_weight) / mb
            return dv_dt, (F_thrust)
            
def equation_vel_air(v, Vw, p_in, _p_ino):
            
            F_drag = Drag(v)
            F_weight = Weight(0)
            T_air = T_initial * (p_in / _p_ino) ** ((gamma - 1) / gamma)
            
            if p_in > patm:
                pressure_ratio = p_in / patm
                
                # Vitesse d'éjection
                v_e = math.sqrt(2 * gamma / (gamma - 1) * R_air * T_air * (1 - (1 / pressure_ratio) ** ((gamma - 1) / gamma)))
                F_thrust = Ae * v_e**2 
                
                
            else:
                F_thrust = 0
            
            dv_dt = (0 - F_drag - F_weight) / mb
            
            return dv_dt, F_thrust

def equation_water_volume(_p_in):

    v_e = water_exit_velocity(_p_in)
    dVw_dt = -Ae * v_e        # rate of change of water volume
    return dVw_dt

#=======================================================
########           Graphique txt             ##########
#=======================================================

# fonction d'affichage des de liste de données suivant l'axe vertical parce que c est plus facil et que ca marche tout aussi bien


def txtGraph(xs: list, ys: list, FileName = "Output.txt"):
    source_file = open(FileName, "w")
    x_simple = xs[:] ##[1, 2.1, 3.6, 4.4, 5.5, 6, 7, 8, 9, 10]  test data
    y_simple = ys[:] ##[2.1, 4, 6, 8, 10, 8.5, 6.6, 4.8, 2.8, 0] test data

    for i in range(len(x_simple)):       # only keep the rounded value of our data because we cant use decimal number to set the position of the text
        x_simple[i] = round(x_simple[i])
        y_simple[i] = round(y_simple[i])

    max_x = max(x_simple)                
    max_y = max(y_simple)                # find the max value of the list to adjust the pos of the text    
    

    for i in range(len(y_simple)):
        line = (y_simple[i] - 1) * " " + str(round(ys[i], 1)) + (max_y - y_simple[i]) * " "
        source_file.write(line + "\n")

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
        Vw = [row[2]*1000 for row in data]
        plt.plot(x, Vw, label="Vw", color="green", marker="o", markersize = 2)

    if len(data[0]) > 3:
        p_in = [row[3]/1000 for row in data]
        plt.plot(x, p_in, label="p_in", color="purple", marker="o", markersize = 2)
    
    if len(data[0]) > 4:
        F = [row[4] for row in data]
        plt.plot(x, F, label="F", color="orange", marker="o", markersize = 2)
    

    
    
    
    
    if additionalValues:
        plt.plot(x, additionalValues, label="F_thrust", color="orange", marker="o", markersize = 2)

    plt.xlabel("Index")
    plt.ylabel("Valeur")
    plt.title("Évolution de h, v et Vw")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()



#=======================================================
########           Runguert-Kutta             ##########
#=======================================================
##### Simulation parameters #####



def systeme_complet(t, y, _p_ino, deltaT, Vw0):
        """
        Params:
        t = time
        y = [h, v, Vw, p_in, F]
        return [dh_dt, dv_dt, dVw_dt, dp_in_dt, F]
        """
        h, v, Vw, p_in, F = y
        _p_ino = _p_ino + patm
        # Équation de position
        dh_dt = v
        p_in_next = internal_pressure(_p_ino, Vw, Vw0)
        # Équation de volume d'eau
        if Vw > 0.0000000001:  # Éviter les valeurs négatives
            dVw_dt = equation_water_volume(p_in)
            
        else:
            dVw_dt = 0
            Vw = 0
        
        # Équation de vitesse
        if Vw > 0.0000000001:
            
            dp_in_dt = (p_in_next - p_in) / deltaT
            dv_dt, F = equation_vel(v, Vw, Vw0, _p_ino)
            
        else:
            if p_in > patm: 
                dp_in_dt = -70000*(1/t)
            else:
                 dp_in_dt = 0

            
            dv_dt, F = equation_vel_air(v, Vw, p_in, _p_ino)

        return [dh_dt, dv_dt, dVw_dt, dp_in_dt, F]

def Rungekutta(_Vw0, _p_ino, stepNbr=50000):
    """
    stepNbr : number of points in the simulation
    """

    # Initial conditions
    t0 = 0.00001 # Initial time
    h0 = 0.00001 # Initial height
    v0 = 0 # Initial velocity
    Vw0 = _Vw0 # Initial water volume
    p_ino = _p_ino # Initial pressure

    mw = rho_w * Vw0
    k = Vw0 /V
    k0 = Vw0/V
    
    
    
    
    # Initialisation
    ts = [t0]
    Fs = [0]
    ys = [[h0, v0, Vw0, p_ino, 0]]  # Stocker toutes les variables dans une liste
    
    deltaT = 0.001  # Pas de temps constant
    i = 0
    
    while ys[-1][0] > 0 and i < stepNbr :
        i += 1
        t_current = ts[-1]
        y_current = ys[-1]
        #print(y_current[0])
        
        p_in = internal_pressure(_p_ino, y_current[2], Vw0)
        v_e = water_exit_velocity(p_in)
        F_thrust = Thrust(v_e)
        Fs.append(p_in)
        #print(p_in)
        # RK4 standard
        k1 = systeme_complet(t_current, y_current, p_ino, deltaT, Vw0)
        
        k2 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k1[j] for j in range(4)]+[y_current[4]], p_ino, deltaT, Vw0)
        
        k3 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k2[j] for j in range(4)]+[y_current[4]], p_ino, deltaT, Vw0)
        
        k4 = systeme_complet(t_current + deltaT, [y_current[j] + deltaT * k3[j] for j in range(4)]+[y_current[4]], p_ino, deltaT, Vw0)
        
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(4)]+[k1[4]]
            
        ts.append(t_current + deltaT)
        ys.append(y_new)


    
    
    #graphMathPlot(ts, ys)
    # Spliting the results
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    Vws = [y[2] for y in ys]
    Fs = [y[4] for y in ys]

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
        print("Veuillez entrer le volume d'eau (en L) :")
        Vw0 = float(input())/1000
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