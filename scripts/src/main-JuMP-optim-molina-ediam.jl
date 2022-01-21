include("utils-JuMP.jl")

#Load parameters required for determining initial conditions
ε = 3
α = 0.33
size_factor = 4
γ_re = 0.25
k_re = 0
γ_ce = 0.25
k_ce = 0
η_re= 0.02
η_ce= 0.02
ν_re = 0.02
ν_ce= 0.02
qsi = 0.010054
δ_S = 0.001823
Δ_T_Disaster= 6
β_T = 4.997053
CO2_base = 289.415046
CO2_Disaster= 1298.216153
labor_growth_N = 0.000
labor_growth_S = 0.000
ρ = 0.015
λ = 0.1443
σ = 2

## Y renewable energy, advanced economies
Yre_N_0 = 45.55074
## Y carbon energy, advanced economies
Yce_N_0 = 193.2
## Y renewable energy, emerging economies
Yre_S_0 = 27.82166
## Y carbon energy, emerging economies
Yce_S_0 = 257.5463
### Environment quality
S_0 = 915.970085

#Initial Productivity conditions are determined by the initial levels of production of energy
#In the Northern Region
global Ace_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_N_0/Yre_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
global Are_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_N_0/Yce_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

#In the Southern Region
global Ace_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_S_0/Yre_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
global Are_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_S_0/Yce_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

# Save de initial state variable in an array
global u0= [Are_N_0, Ace_N_0, Are_S_0,Ace_S_0,S_0]

# Save parameters in an array
global dt = 5    # 1 Quarterly
global D = 120.0 # Simulate for 30 years

model = Model()
model = optimize_model(model,"P2")
solution_summary(model, verbose=true)
