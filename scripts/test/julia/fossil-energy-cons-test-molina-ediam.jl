using Turing
using DifferentialEquations
using Plots,StatsPlots


function molina_ediam(du,u,p,t)
    ## Parámetros iniciales
    ε,α,size_factor,γ_re,k_re,γ_ce,k_ce,η_re,η_ce,ν_re,ν_ce,qsi,δ_S,Δ_T_Disaster,β_T,CO2_base,CO2_Disaster,labor_growth_N,labor_growth_S,ρ,λ,σ = p

    #=
    Derivada con respecto al tiempo del vector de estado.
        * u es el vector de estado (arreglo)
    =#
    time =  t

    Are_N,Ace_N,Are_S,Ace_S,S = u

    Are_N_0=83.95750622106519
    Ace_N_0 =177.31606561865874

    ### Auxiliares generales

    φ= (1-α)*(1-ε)

    #this is the cost of production of clean technologies
    epsi_re = α^2
    #this is the cost of production of dirty technologies
    epsi_ce = α^2

    ### North Region
    #Auxiliaries in North

    L_N = ℯ^(labor_growth_N*time)

    #gamma displays decreasing returns as in Stiligtz
    γ_re_t_N = γ_re*ℯ^(-k_re*(Are_N/Are_N_0-1))

    #gamma displays decreasing returns as in Stiligtz
    γ_ce_t_N = γ_ce*ℯ^(-k_ce*(Ace_N/Ace_N_0-1))

    ### Carbon tax in advanced region
    ce_tax_N=0
    ### Technology subsidy in advanced region
    Tec_subsidy_N=0

    ### Subsidies for research and development
    RD_subsidy_N = 0

    RelPrice_N = ((Ace_N/Are_N)^(1-α))*(((epsi_re*(1-Tec_subsidy_N))/epsi_ce)^α)
    RelLabor_N =((1+ce_tax_N)^ε)*((((1-Tec_subsidy_N)*epsi_re)/epsi_ce)^(α*(1-ε)))*((Are_N/Ace_N)^(-1*φ))

    # Clean sector
    #based on the assumption that Labor.re.N+Labor.ce.N=L.N
    Labor_re_N = (RelLabor_N*L_N)/(1+RelLabor_N)
    #based on the assumption that  Price.re.N^(1-ε)+Price.ce.N^(1-ε)=1
    Price_re_N = RelPrice_N/(RelPrice_N^(1-ε)+(1)^(1-ε))^(1/(1-ε))
    # technology demand
    Agg_demand_re_tech_N = ((((α^2)*Price_re_N)/((1-Tec_subsidy_N)*epsi_re))^(1/(1-α)))*Labor_re_N*Are_N
    # Expected profits see annex IV. Equilibrium research profits
    Profits_re_N =(1+RD_subsidy_N)*η_re*epsi_re*((1-α)/α)*Agg_demand_re_tech_N
    # Equilibrium levels of production
    Yre_N = ((((α^2)*Price_re_N)/((1-Tec_subsidy_N)*epsi_re))^(α/(1-α)))*Labor_re_N*Are_N

    # dirty sector
    Labor_ce_N = L_N/(RelLabor_N+1)
    Price_ce_N = Price_re_N/RelPrice_N
    Agg_demand_ce_tech_N = ((((α^2)*Price_ce_N)/(epsi_ce))^(1/(1-α)))*Labor_ce_N*Ace_N
    Profits_ce_N = η_ce*epsi_ce*((1-α)/α)*Agg_demand_ce_tech_N
    Yce_N = ((((α^2)*Price_ce_N)/(epsi_ce))^(α/(1-α)))*Labor_ce_N*Ace_N

    # Producción total

    Y_N = ((Yre_N)^((ε-1)/ε)+(Yce_N)^((ε-1)/ε))^(ε/(ε-1))

    sre_N = ℯ^(Profits_re_N)/(ℯ^(Profits_ce_N)+ℯ^(Profits_re_N))
    sce_N = 1-sre_N

    #Auxiliaries in South
    #the population of the South is 4.6 that of the North,
    L_S = (ℯ^(labor_growth_S*time))*size_factor
    γ_re_t_S = γ_re
    γ_ce_t_S = γ_ce

    ### Carbon tax in emergent region
    ce_tax_S=0
    ### Technology subsidy in emergent region
    Tec_subsidy_S=0

    ### Subsidies for research and development
    RD_subsidy_S = 0
    #First we determine the equilibrium levels of relative input prices and relative labour

    RelPrice_S = ((Ace_S/Are_S)^(1-α))*(((epsi_re*(1-Tec_subsidy_S))/epsi_ce)^α)
    RelLabor_S = ((1+ce_tax_S)^ε)*((((1-Tec_subsidy_S)*epsi_re)/epsi_ce)^(α*(1-ε)))*((Are_S/Ace_S)^(-1*φ))

    #Second we determine the equilibrium conditions for each sector
    #clean sector
    #based on the assumption that Labor_re_S+Labor_ce_S=L_S
    Labor_re_S = (L_S*RelLabor_S)/(RelLabor_S+1)
    #based on the assumption that  Price_re_S^(1-ε)+(Price_ce_S)^(1-ε)=1
    Price_re_S = RelPrice_S/(RelPrice_S^(1-ε)+(1)^(1-ε))^(1/(1-ε))
    Agg_demand_re_tech_S = ((((α^2)*Price_re_S)/((1-Tec_subsidy_S)*epsi_re))^(1/(1-α)))*Labor_re_S*Are_S
    Profits_re_S = (1+RD_subsidy_S)*η_re*epsi_re*((1-α)/α)*Agg_demand_re_tech_S
    Yre_S = ((((α^2)*Price_re_S)/((1-Tec_subsidy_S)*epsi_re))^(α/(1-α)))*Labor_re_S*Are_S

    #dirty sector
    Labor_ce_S = L_S/(RelLabor_S+1)
    Price_ce_S = Price_re_S/RelPrice_S
    Agg_demand_ce_tech_S = ((((α^2)*Price_ce_S)/(epsi_ce))^(1/(1-α)))*Labor_ce_S*Ace_S
    Profits_ce_S = η_ce*epsi_ce*((1-α)/α)*Agg_demand_ce_tech_S
    Yce_S = ((((α^2)*Price_ce_S)/(epsi_ce))^(α/(1-α)))*Labor_ce_S*Ace_S

    #Total Production
    Y_S = ((Yre_S)^((ε-1)/ε)+(Yce_S)^((ε-1)/ε))^(ε/(ε-1))

    #Allocation of Scientists
    sre_S = ℯ^(Profits_re_S)/(ℯ^(Profits_ce_S)+ℯ^(Profits_re_S))
    sce_S = 1-sre_S

    ##### Changes in Temperature
    #increase in temperature at which there is environmental disaster
    Delta_Temp_Disaster = Δ_T_Disaster
    CO2_Concentration = max(CO2_Disaster-S,CO2_base)
    Delta_Temp = min(β_T*log(CO2_Concentration/CO2_base),Delta_Temp_Disaster)

    #Welfare Calculations
    Consumption_N = Y_N-epsi_re*Agg_demand_re_tech_N-epsi_ce*Agg_demand_ce_tech_N
    Consumption_S = (Y_S-epsi_re*Agg_demand_re_tech_S-epsi_ce*Agg_demand_ce_tech_S)*(1/size_factor)
    Cost_S_Damage = ((Delta_Temp_Disaster-Delta_Temp)^λ-λ*Delta_Temp_Disaster^(λ-1)*(Delta_Temp_Disaster-Delta_Temp))/((1-λ)*Delta_Temp_Disaster^λ)

    #Budget restrictions
    Tec_subsidy_GF_N = 0
    RD_subsidy_GF_N = 0
    Budget_function_N = (ce_tax_N*Price_ce_N*Yce_N) - (Tec_subsidy_N*epsi_re*Agg_demand_re_tech_N) - (Tec_subsidy_GF_N*epsi_re*Agg_demand_re_tech_S) -(RD_subsidy_N*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_N )-(RD_subsidy_GF_N*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_S)

    Budget_function_S = (ce_tax_S*Price_ce_S*Yce_S)- (Tec_subsidy_S*epsi_re*Agg_demand_re_tech_S) - (RD_subsidy_S*η_re*((epsi_re/α)-epsi_re)*Agg_demand_re_tech_S)

    # Save non-state variables
    push!(N_fossil_energy,Yce_N)
    push!(S_fossil_energy,Yce_S)
    push!(r_Cost_S_Damage,Cost_S_Damage)
    push!(r_Consumption_N,Consumption_N)
    push!(r_Consumption_S,Consumption_S)
    push!(r_Delta_Temp,Delta_Temp)

    #State variables
    #Evolution of Productivity North Region
    dAre_N = γ_re_t_N*η_re*sre_N*Are_N
    dAce_N = γ_ce_t_N*η_ce*sce_N*Ace_N

    #Evolution of Productivity South Region
    dAre_S = γ_re_t_S*ν_re*sre_S*(Are_N-Are_S)
    dAce_S = γ_ce_t_S*ν_ce*sce_S*(Ace_N-Ace_S)

    #Environmental Quality
    dS = min(1_0,δ_S*S-qsi*(Yce_N+Yce_S))

    du .=(dAre_N,dAce_N,dAre_S,dAce_S,dS)
end


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


p = [ε ,α ,size_factor,γ_re ,k_re,γ_ce ,k_ce,η_re,η_ce,ν_re ,ν_ce,qsi ,δ_S ,Δ_T_Disaster,β_T ,CO2_base ,CO2_Disaster,labor_growth_N ,labor_growth_S,ρ,λ,σ]

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
Ace_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_N_0/Yre_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
Are_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_N_0/Yce_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

#In the Southern Region
Ace_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_S_0/Yre_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
Are_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_S_0/Yce_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))


u0= [Are_N_0, Ace_N_0, Are_S_0,Ace_S_0,S_0]

# Define array for the non-state variables
N_fossil_energy = []
S_fossil_energy = []
r_Cost_S_Damage = []
r_Consumption_N = []
r_Consumption_S = []
r_Delta_Temp = []

prob1 = ODEProblem(molina_ediam,u0,(1.0,30.0),p)
sol = solve(prob1,Euler(),dt=0.25)

# Graficamos el Consumo de energia fósil en ambas regiones
plot([1983+i for i in sol.t],hcat(N_fossil_energy,S_fossil_energy),
 title = "Consumo de energia fósil",
 label = ["Regiones Avanzadas" "Regiones Emergentes"],
  lw = 3)

# Use simulation output to estimate value of objective function
# (1+((out$Cost.S.Damage*out$Consumption_N)^(1-sigma.utility))/(1-sigma.utility))*(1/((1+rho)^out$time))
utility_consumer_N = 1 .+ ((r_Cost_S_Damage .* r_Consumption_N).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))
utility_consumer_S = 1 .+ ((r_Cost_S_Damage .* r_Consumption_S).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))

# Greficamos el incremento de la temperatura
plot([1983+i for i in sol.t],r_Delta_Temp,
 title = "Incremento de la temperatura",
  lw = 3)


# https://discourse.julialang.org/t/using-differentialequations-jl-to-solve-routine-optimal-control-problems/34096
#https://discourse.julialang.org/t/solving-odes-in-a-loop-with-differentialequations-jl/59640
#https://nextjournal.com/sosiris-de/ode-diffeq?change-id=CkQATVFdWBPaEkpdm6vuto
#https://discourse.julialang.org/t/differentialequations-jl-discrete-delay-problem/62829/11
#https://diffeq.sciml.ai/stable/types/discrete_types/
#https://sciml.ai/news/2018/04/30/Jupyter/
#https://discourse.julialang.org/t/differentialequations-inexact-error/31208
#https://discourse.julialang.org/t/solving-difference-equation/66977/2
#https://discourse.julialang.org/t/solving-difference-equation/66977/3
#https://github.com/mitmath/18337
# http://rstudio-pubs-static.s3.amazonaws.com/32888_197d1a1896534397b67fb04e0d4899ae.html