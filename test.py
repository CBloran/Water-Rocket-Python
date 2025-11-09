import math


# ============================================================
# CONSTANTES PHYSIQUES RÉALISTES
# ============================================================

g = 9.81               # Gravité (m/s²)
rho_w = 1000.0         # Densité eau (kg/m³)
rho_air = 1.225        # Densité air (kg/m³)
Cd = 0.35               # Coefficient traînée (fusée typique)
patm = 101325          # Pression atmosphérique (Pa)

# Paramètres FUSÉE RÉALISTE
D_body = 0.1           # Diamètre fusée (m)
A_cross = math.pi * (D_body/2)**2  # Section frontale
D_nozzle = 0.02       # Diamètre buse (m) - PLUS PETIT
A_nozzle = math.pi * (D_nozzle/2)**2

V_total = 0.0015        # Volume total 2L
V_water_initial = 0.0005 # Volume eau initial 1L
V_air_initial = V_total - V_water_initial

p0 = 500000           # Pression initiale 4 bars
m_dry = 0.1           # Masse à vide 500g - PLUS LÉGER

print("=== PARAMÈTRES FUSÉE ===")
print(f"Volume eau: {V_water_initial*1000:.0f}mL")
print(f"Volume air: {V_air_initial*1000:.0f}mL") 
print(f"Pression: {p0/1000:.0f} kPa")
print(f"Masse sèche: {m_dry:.1f}kg")
print(f"Surface buse: {A_nozzle*10000:.1f} cm²")

class WaterRocket:
    def __init__(self):
        self.g = g
        self.rho_w = rho_w
        self.rho_air = rho_air
        self.Cd = Cd
        self.patm = patm
        self.A_cross = A_cross
        self.A_nozzle = A_nozzle
        self.V_total = V_total
        self.V_water = V_water_initial
        self.V_air = V_air_initial
        self.pressure = p0
        self.m_dry = m_dry
        self.m_water = rho_w * V_water_initial
        
    def update(self, dt):
        """Mise à jour physique sur un pas de temps dt"""
        
        # 1. ÉJECTION DE L'EAU (si il reste de l'eau et pression suffisante)
        thrust = 0
        if self.V_water > 0 and self.pressure > self.patm:
            # Vitesse d'éjection (Bernoulli)
            v_ejection = math.sqrt(2 * (self.pressure - self.patm) / self.rho_w)
            
            # Débit volumique
            Q = self.A_nozzle * v_ejection
            
            # Masse d'eau éjectée
            dm_water = self.rho_w * Q * dt
            dm_water = min(dm_water, self.m_water)  # Ne pas éjecter plus que disponible
            
            if dm_water > 0:
                # Mise à jour masse et volume eau
                self.m_water -= dm_water
                self.V_water = self.m_water / self.rho_w
                
                # Force de poussée (F = dm/dt * v)
                thrust = (dm_water / dt) * v_ejection
                
                # Expansion adiabatique de l'air
                self.V_air = self.V_total - self.V_water
                if self.V_air > 0:
                    # Loi adiabatique: P * V^γ = constante
                    gamma = 1.4
                    self.pressure = p0 * (V_air_initial / self.V_air) ** gamma
                else:
                    self.pressure = self.patm
        
        return thrust
    
    def get_mass(self):
        return self.m_dry + self.m_water
    
    def get_drag(self, velocity):
        """Force de traînée"""
        if abs(velocity) < 0.1:
            return 0
        return 0.5 * self.rho_air * velocity**2 * self.Cd * self.A_cross

def simulate_rocket():
    """Simulation complète de la trajectoire"""
    
    rocket = WaterRocket()
    
    # Conditions initiales
    t = 0
    dt = 0.01  # Pas de temps 10ms
    y = 0.1    # Hauteur initiale (départ du sol)
    v = 0.0    # Vitesse initiale
    phase = "PROPULSION"
    
    # Stockage résultats
    times = [t]
    heights = [y]
    velocities = [v]
    thrusts = [0]
    masses = [rocket.get_mass()]
    pressures = [rocket.pressure]
    
    max_height = 0
    water_depleted_time = 0
    
    print("\n=== DÉBUT SIMULATION ===")
    
    while y > 0 or t < 1:  # Simuler jusqu'au sol ou minimum 1s
        # Calcul forces
        thrust = rocket.update(dt)
        drag = rocket.get_drag(v)
        weight = rocket.get_mass() * g
        
        # Accélération (F = ma)
        if rocket.get_mass() > 0.01:  # Éviter division par zéro
            acceleration = (thrust - drag - weight) / rocket.get_mass()
        else:
            acceleration = -g  # Chute libre
            
        # Intégration vitesse et position
        v += acceleration * dt
        y += v * dt
        
        # Mise à jour temps
        t += dt
        
        # Détection fin eau
        if rocket.V_water <= 0.000001 and water_depleted_time == 0:
            water_depleted_time = t
            phase = "BALLISTIQUE"
            print(f"→ Phase balistique à t={t:.2f}s, h={y:.1f}m, v={v:.1f}m/s")
        
        # Stockage données
        times.append(t)
        heights.append(max(y, 0))
        velocities.append(v)
        thrusts.append(thrust)
        masses.append(rocket.get_mass())
        pressures.append(rocket.pressure)
        
        max_height = max(max_height, y)
        
        # Arrêt si au sol depuis un moment
        if y <= 0 and t > 2:
            break
            
        # Sécurité durée
        if t > 30:
            break
    
    # Analyse résultats
    print("\n=== RÉSULTATS ===")
    print(f"Temps simulation: {t:.2f}s")
    print(f"Hauteur maximale: {max_height:.1f}m")
    print(f"Vitesse maximale: {max([abs(v) for v in velocities]):.1f}m/s")
    print(f"Poussée max: {max(thrusts):.1f}N")
    print(f"Temps propulsion: {water_depleted_time:.2f}s")
    
    if max_height < 5:
        print("\n⚠️  HAUTEUR TROP FAIBLE - VÉRIFIER PARAMÈTRES")
        print("Suggestions:")
        print("- Augmenter pression initiale")
        print("- Réduire masse sèche")
        print("- Augmenter volume eau")
        print("- Agrandir diamètre buse")
    
    return times, heights, velocities, thrusts, masses, pressures

def create_ascii_trajectory(times, heights):
    """Crée une visualisation ASCII de la trajectoire"""
    if not heights:
        return
    
    max_h = max(heights)
    scale_h = 50.0 / max_h if max_h > 0 else 1
    scale_t = len(times) / 50  # Échantillonnage temporel
    
    print(f"\n📈 TRAJECTOIRE (max: {max_h:.1f}m)")
    print("=" * 60)
    
    for i in range(0, len(times), max(1, int(scale_t))):
        if i < len(heights) and heights[i] >= 0:
            h_display = heights[i] * scale_h
            bar = "█" * int(h_display)
            print(f"t={times[i]:5.2f}s | h={heights[i]:5.1f}m {bar}")

# ============================================================
# EXÉCUTION
# ============================================================

if __name__ == "__main__":
    print("SIMULATION FUSÉE À EAU")
    print("Modèle physique réaliste avec propulsion par éjection d'eau")
    
    times, heights, velocities, thrusts, masses, pressures = simulate_rocket()
    
    create_ascii_trajectory(times, heights)
    
    # Sauvegarde données
    with open("rocket_trajectory.txt", "w") as f:
        f.write("t(s)\th(m)\tv(m/s)\tF(N)\tm(kg)\tp(Pa)\n")
        for i in range(len(times)):
            f.write(f"{times[i]:.3f}\t{heights[i]:.3f}\t{velocities[i]:.3f}\t")
            f.write(f"{thrusts[i]:.1f}\t{masses[i]:.3f}\t{pressures[i]:.0f}\n")
    
    print(f"\n💾 Données sauvegardées dans 'rocket_trajectory.txt'")