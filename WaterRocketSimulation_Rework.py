
import math
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

mb = 0.058
rho_eau = 1000
P_atm = 101325
Cd = 0.71
rho_air = 1.225
D_tuyere = 0.009
D_fusee = 0.09
A_tuyere = math.pi * (D_tuyere**2) / 4
A_frontale =  math.pi * (D_fusee**2) / 4
V = 0.0015
g = 9.81





def vitesseEjection(P):
 
    delta_P = max(P - P_atm, 0.0)
 
    return -math.sqrt(2.0 * delta_P / (rho_eau*(1-(A_tuyere/A_frontale)**2)))
 

def calculateF_grav(mw):
    F_grav = (mb + mw) * g
    if mw > 0.001:
        print(F_grav)
    return F_grav

def calculateDrag(v):

    return 0.5 * rho_air * Cd * A_frontale * v **2

def Pression(mw0, p0, mw):
    
        #P = (p0*((V - mw0/rho_eau)**1.4))/((V - mw/rho_eau)**1.4)
        P = p0 * ((V - mw0/rho_eau)**1.4)/((V - mw/rho_eau)**1.4)
        return P




def Euler(mw0, p0):

    
    
    h0 = 0
    v0 = 0
    ts = [0]

    ys = [[h0, v0, mw0, 0]]
    ts = [0]
    deltaT = 0.001

    
    def systemeCompletEau(y):
        
        h, v, mw, = y

        dh_dt = v
        
        P = Pression(mw0, p0, mw)
        F_grav = calculateF_grav(mw)
        F_trainee = calculateDrag(v)
        v_e = vitesseEjection(P)

        dmw_dt = rho_eau * A_tuyere * v_e

        F_poussee = dmw_dt*v_e

        dv_dt = (F_poussee - F_grav - F_trainee)/(mb + mw) 

        

        return [dh_dt, dv_dt, dmw_dt], (F_poussee - F_grav)
    
    def systemeCompletAutre(y):
        
        h, v, mw, = y

        dh_dt = v
        
        P = P_atm
        F_grav = calculateF_grav(mw)
        F_trainee = calculateDrag(v)
        v_e = 0

        F_poussee = 0

        dv_dt = (F_poussee - F_grav - F_trainee)/(mb) 

        dmw_dt = 0

        return [dh_dt, dv_dt, dmw_dt], (-F_grav)


    while ys[-1][0] >= 0:
        
        t_current = ts[-1]
        y_current = ys[-1]
        
        CalculatedVal = []

        if y_current[2] > 0.00000001:
                
            k, F = systemeCompletEau(y_current[:3])
            CalculatedVal.append(F)
            y_new = [y_current[i] + deltaT * k[i] for i in range(3)] + CalculatedVal
            #print(y_current)
            ys.append(y_new)
            ts.append(t_current + deltaT)
        
        else:
            k, F = systemeCompletAutre(y_current[:3])
            CalculatedVal.append(F)
            y_new = [y_current[i] + deltaT * k[i] for i in range(3)] + CalculatedVal
            
            ys.append(y_new)
            ts.append(t_current + deltaT)
    
    with open("resultats.txt", "w", encoding="utf-8")as w:
        for i in range(len(ys)):
            line =  str(ts[i])
            for val in ys[i]:
                line += "," + str(val)
            line += "\n"
            w.write(line)
    return ts, ys
    

def EulerMax(x, y):
    result = Euler(x, y)[1]
    hs = [y[0] for y in result]
    maxh = max(hs)
    
    return maxh


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
        Vw = [row[2] for row in data]
        plt.plot(x, Vw, label="mw", color="green", marker="o", markersize = 2)

    if len(data[0]) > 3:
        F = [row[3] for row in data]
        plt.plot(x, F, label="F", color="purple", marker="o", markersize = 2)
    
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


def graphique_3d_parametres(fonction, x_min, x_max, y_min, y_max, nb_points_x=40, nb_points_y=40):
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
    ax.set_xlabel("Masse d'eau (Kg))")
    ax.set_ylabel("Pression d'entrée (Pa)")
    ax.set_zlabel('Hauteur atteinte')
    ax.set_title("Evolution de la hauteur atteinte en fonction du volume d'eau et de la pression")
    
    # Ajout d'une barre de couleur
    plt.colorbar(surface, ax=ax, shrink=0.5, aspect=5)
    
    plt.show()
    
    return fig, ax
#graphique_3d_parametres(RungekuttaMax, 0.0001, 0.0010, 100000, 500000)

def batch(mw0):
    array = []
    i = 100000
    while i <= 500000:
        print(i)
        hmax = EulerMax(mw0, i + P_atm)
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
        print("Veuillez entrer le masse d'eau (en Kg) :")
        mw0 = float(input())
        print("Veuillez entrer la pression d'entrée (en bar) :")
        p_in0 = float(input())*10**5 + P_atm
        
        
        tsys = Euler(mw0, p_in0)
        print("Hauteur maximale atteinte : ", max(tsys[1]))
        graphMathPlot(tsys[0], tsys[1])
    elif choice == 2:
        print("!!WARNING THIS CAN TAKE A WHILE !!")
        print("Press any key to continue")
        input()
        graphique_3d_parametres(EulerMax, 0.1, 1, 1000, 5*10**5)
    elif choice == 3:
        mw0 = float(input("Veuillez entrer le masse d'eau (en Kg) : "))
        batch(mw0)
    elif choice == 4:
        print("Veuillez entrer la hauteur maximale souhaitée (en m) :")
        Height = float(input())
        print("La pression d'entrée optimale pour atteindre une hauteur de ", Height, "m est de : ", FinBestPressure(Height), "Pa")
        