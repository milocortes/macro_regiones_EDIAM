import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

def ode_FE(f, U_0, dt, T):
    N_t = int(round(T/dt))
    # Ensure that any list/tuple returned from f_ is wrapped as array
    f_ = lambda u, t: np.asarray(f(u, t))
    u = np.zeros((N_t+1, len(U_0)))
    t = np.linspace(1, N_t*dt, len(u))

    u[0] = U_0

    for n in range(N_t):
        u[n+1] = u[n] + dt*f_(u[n], t[n])

    return u, t


def f(u, t):

    time = t

    Are_N,Ace_N,Are_S,Ace_S,S = u

    ### Auxiliares generales

    φ= (1-α)*(1-ε)

    #this is the cost of production of clean technologies
    epsi_re = α**2
    #this is the cost of production of dirty technologies
    epsi_ce = α**2

    ### North Region
    #Auxiliaries in North

    L_N = math.exp(labor_growth_N*time)

    #gamma displays decreasing returns as in Stiligtz
    γ_re_t_N = γ_re*math.exp(-k_re*(Are_N/Are_N_0-1))

    #gamma displays decreasing returns as in Stiligtz
    γ_ce_t_N = γ_ce*math.exp(-k_ce*(Ace_N/Ace_N_0-1))

    ### Carbon tax in advanced region
    ce_tax_N=0
    ### Technology subsidy in advanced region
    Tec_subsidy_N=0

    ### Subsidies for research and development
    RD_subsidy_N = 0

    RelPrice_N = ((Ace_N/Are_N)**(1-α))*(((epsi_re*(1-Tec_subsidy_N))/epsi_ce)**α)
    RelLabor_N =((1+ce_tax_N)**ε)*((((1-Tec_subsidy_N)*epsi_re)/epsi_ce)**(α*(1-ε)))*((Are_N/Ace_N)**(-1*φ))

    # Clean sector
    #based on the assumption that Labor.re.N+Labor.ce.N=L.N
    Labor_re_N = (RelLabor_N*L_N)/(1+RelLabor_N)
    #based on the assumption that  Price.re.N**(1-ε)+Price.ce.N**(1-ε)=1
    Price_re_N = RelPrice_N/(RelPrice_N**(1-ε)+(1)**(1-ε))**(1/(1-ε))
    # technology demand
    Agg_demand_re_tech_N = ((((α**2)*Price_re_N)/((1-Tec_subsidy_N)*epsi_re))**(1/(1-α)))*Labor_re_N*Are_N
    # Expected profits see annex IV. Equilibrium research profits
    Profits_re_N =(1+RD_subsidy_N)*η_re*epsi_re*((1-α)/α)*Agg_demand_re_tech_N
    # Equilibrium levels of production
    Yre_N = ((((α**2)*Price_re_N)/((1-Tec_subsidy_N)*epsi_re))**(α/(1-α)))*Labor_re_N*Are_N

    # dirty sector
    Labor_ce_N = L_N/(RelLabor_N+1)
    Price_ce_N = Price_re_N/RelPrice_N
    Agg_demand_ce_tech_N = ((((α**2)*Price_ce_N)/(epsi_ce))**(1/(1-α)))*Labor_ce_N*Ace_N
    Profits_ce_N = η_ce*epsi_ce*((1-α)/α)*Agg_demand_ce_tech_N
    Yce_N = ((((α**2)*Price_ce_N)/(epsi_ce))**(α/(1-α)))*Labor_ce_N*Ace_N

    # Producción total

    Y_N = ((Yre_N)**((ε-1)/ε)+(Yce_N)**((ε-1)/ε))**(ε/(ε-1))

    sre_N = math.exp(Profits_re_N)/(math.exp(Profits_ce_N)+math.exp(Profits_re_N))
    sce_N = 1-sre_N

    #Auxiliaries in South
    #the population of the South is 4.6 that of the North,
    L_S = (math.exp(labor_growth_S*time))*size_factor
    γ_re_t_S = γ_re
    γ_ce_t_S = γ_ce

    ### Carbon tax in emergent region
    ce_tax_S=0
    ### Technology subsidy in emergent region
    Tec_subsidy_S=0

    ### Subsidies for research and development
    RD_subsidy_S = 0
    #First we determine the equilibrium levels of relative input prices and relative labour

    RelPrice_S = ((Ace_S/Are_S)**(1-α))*(((epsi_re*(1-Tec_subsidy_S))/epsi_ce)**α)
    RelLabor_S = ((1+ce_tax_S)**ε)*((((1-Tec_subsidy_S)*epsi_re)/epsi_ce)**(α*(1-ε)))*((Are_S/Ace_S)**(-1*φ))

    #Second we determine the equilibrium conditions for each sector
    #clean sector
    #based on the assumption that Labor_re_S+Labor_ce_S=L_S
    Labor_re_S = (L_S*RelLabor_S)/(RelLabor_S+1)
    #based on the assumption that  Price_re_S**(1-ε)+(Price_ce_S)**(1-ε)=1
    Price_re_S = RelPrice_S/(RelPrice_S**(1-ε)+(1)**(1-ε))**(1/(1-ε))
    Agg_demand_re_tech_S = ((((α**2)*Price_re_S)/((1-Tec_subsidy_S)*epsi_re))**(1/(1-α)))*Labor_re_S*Are_S
    Profits_re_S = (1+RD_subsidy_S)*η_re*epsi_re*((1-α)/α)*Agg_demand_re_tech_S
    Yre_S = ((((α**2)*Price_re_S)/((1-Tec_subsidy_S)*epsi_re))**(α/(1-α)))*Labor_re_S*Are_S

    #dirty sector
    Labor_ce_S = L_S/(RelLabor_S+1)
    Price_ce_S = Price_re_S/RelPrice_S
    Agg_demand_ce_tech_S = ((((α**2)*Price_ce_S)/(epsi_ce))**(1/(1-α)))*Labor_ce_S*Ace_S
    Profits_ce_S = η_ce*epsi_ce*((1-α)/α)*Agg_demand_ce_tech_S
    Yce_S = ((((α**2)*Price_ce_S)/(epsi_ce))**(α/(1-α)))*Labor_ce_S*Ace_S

    #Total Production
    Y_S = ((Yre_S)**((ε-1)/ε)+(Yce_S)**((ε-1)/ε))**(ε/(ε-1))

    #Allocation of Scientists
    sre_S = math.exp(Profits_re_S)/(math.exp(Profits_ce_S)+math.exp(Profits_re_S))
    sce_S = 1-sre_S

    ##### Changes in Temperature
    #increase in temperature at which there is environmental disaster
    Delta_Temp_Disaster = Δ_T_Disaster
    CO2_Concentration = max(CO2_Disaster-S,CO2_base)
    Delta_Temp = min(β_T*math.log(CO2_Concentration/CO2_base),Delta_Temp_Disaster)
    Delta_Temp_list.append(Delta_Temp)

    #Welfare Calculations
    Consumption_N = Y_N-epsi_re*Agg_demand_re_tech_N-epsi_ce*Agg_demand_ce_tech_N
    Consumption_S = (Y_S-epsi_re*Agg_demand_re_tech_S-epsi_ce*Agg_demand_ce_tech_S)*(1/size_factor)
    Cost_S_Damage = ((Delta_Temp_Disaster-Delta_Temp)**λ-λ*Delta_Temp_Disaster**(λ-1)*(Delta_Temp_Disaster-Delta_Temp))/((1-λ)*Delta_Temp_Disaster**λ)

    #Budget restrictions
    Tec_subsidy_GF_N = 0
    RD_subsidy_GF_N = 0
    Budget_function_N = (ce_tax_N*Price_ce_N*Yce_N) - (Tec_subsidy_N*epsi_re*Agg_demand_re_tech_N) - (Tec_subsidy_GF_N*epsi_re*Agg_demand_re_tech_S) -(RD_subsidy_N*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_N )-(RD_subsidy_GF_N*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_S)

    Budget_function_S = (ce_tax_S*Price_ce_S*Yce_S)- (Tec_subsidy_S*epsi_re*Agg_demand_re_tech_S) - (RD_subsidy_S*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_S)

    N_fossil_energy.append(Yce_N)
    S_fossil_energy.append(Yce_S)

    return [γ_re_t_N*η_re*sre_N*Are_N,γ_ce_t_N*η_ce*sce_N*Ace_N,γ_re_t_S*ν_re*sre_S*(Are_N-Are_S),γ_ce_t_S*ν_ce*sce_S*(Ace_N-Ace_S),min(1_0,δ_S*S-qsi*(Yce_N+Yce_S))]

### Molina EDIAM
#Load parameters required for determining initial conditions
ε = 3.5
α = 0.33
size_factor = 1
γ_re = 0.25
k_re = 0
γ_ce = 0.25
k_ce = 0
η_re= 0.03
η_ce= 0.021
ν_re = 0.3
ν_ce= 0.25
qsi = 0.010054
δ_S = 0.001823
Δ_T_Disaster= 7.5
β_T = 4.997053
CO2_base = 289.415046
CO2_Disaster= 1298.216153
labor_growth_N = 0.008
labor_growth_S = 0.015
ρ = 0.01
λ = 0.1443
σ = 2
ξ = 0.01

## Y renewable energy, advanced economies
Yre_N_0 = 25.1
## Y carbon energy, advanced economies
Yce_N_0 = 144.9
## Y renewable energy, emerging economies
Yre_S_0 = 9.0
## Y carbon energy, emerging economies
Yce_S_0 = 105.3
### Environment quality
S_0 = 915.970085

#Initial Productivity conditions are determined by the initial levels of production of energy
#In the Northern Region
Ace_N_0 = ((Yce_N_0**((ε-1)/ε)+Yre_N_0**((ε-1)/ε))**(ε/(ε-1)))*(1+(Yce_N_0/Yre_N_0)**((1-ε)/ε))**(1/((1-α)*(1-ε)))
Are_N_0 = ((Yce_N_0**((ε-1)/ε)+Yre_N_0**((ε-1)/ε))**(ε/(ε-1)))*(1+(Yre_N_0/Yce_N_0)**((1-ε)/ε))**(1/((1-α)*(1-ε)))

#In the Southern Region
Ace_S_0 = (1/size_factor)*((Yce_S_0**((ε-1)/ε)+Yre_S_0**((ε-1)/ε))**(ε/(ε-1)))*(1+(Yce_S_0/Yre_S_0)**((1-ε)/ε))**(1/((1-α)*(1-ε)))
Are_S_0 = (1/size_factor)*((Yce_S_0**((ε-1)/ε)+Yre_S_0**((ε-1)/ε))**(ε/(ε-1)))*(1+(Yre_S_0/Yce_S_0)**((1-ε)/ε))**(1/((1-α)*(1-ε)))

dt = 0.1    # 1 Quarterly
D = 30 # Simulate for 30 years
N_t = int(D/dt) # Corresponding no of time steps
T = dt*N_t    # End time
U_0 = [Are_N_0, Ace_N_0, Are_S_0,Ace_S_0,S_0]

N_fossil_energy = []
S_fossil_energy = []
Delta_Temp_list = []
u, t = ode_FE(f, U_0, dt, T)
#N_fossil_energy = u[:,1]
#S_fossil_energy = u[:,3]
fig = plt.figure()
tiempo_line = [1983 +x for x in t[:-1]]
l1, l2 = plt.plot(tiempo_line, N_fossil_energy, tiempo_line, S_fossil_energy)
fig.legend((l1, l2), ("Regiones Avanzadas", "Regiones Emergentes"), "center right")
plt.xlabel("Meses")
plt.show()

if __name__ == '__main__':
    demo_molina()
