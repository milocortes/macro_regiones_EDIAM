using Turing
using DifferentialEquations
using Plots,StatsPlots
using Optim

function molina_ediam(du,u,p,t)
    ## Parámetros iniciales
    ε,α,size_factor,γ_re,k_re,γ_ce,k_ce,η_re,η_ce,ν_re,ν_ce,qsi,δ_S,Δ_T_Disaster,β_T,CO2_base,CO2_Disaster,labor_growth_N,labor_growth_S,ρ,λ,σ,index_vector,ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S = p

    #=
    Derivada con respecto al tiempo del vector de estado.
        * u es el vector de estado (arreglo)
    =#
    time =  t

    index = findall(x -> x == time,index_vector)

    Are_N,Ace_N,Are_S,Ace_S,S = u

    Are_N_0 = 163.31730041748838
    Ace_N_0 = 245.3275202070854

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
    #ce_tax_N = ce_tax_N * ce_tax_N_int[index][1]
    ### Technology subsidy in advanced region
    #Tec_subsidy_N=0

    ### Subsidies for research and development
    #RD_subsidy_N = 0

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
    #ce_tax_S = ce_tax_S * ce_tax_S_int[index][1]
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
    push!(tiempo_line,time)

    #State variables
    #Evolution of Productivity North Region
    dAre_N = γ_re_t_N*η_re*sre_N*Are_N
    dAce_N = γ_ce_t_N*η_ce*sce_N*Ace_N

    #Evolution of Productivity South Region
    dAre_S = γ_re_t_S*ν_re*sre_S*(Are_N-Are_S)
    dAce_S = γ_ce_t_S*ν_ce*sce_S*(Ace_N-Ace_S)

    #Environmental Quality
    dS = min(1.0,(δ_S*S)-qsi*(Yce_N+Yce_S))

    du .=(dAre_N,dAce_N,dAre_S,dAce_S,dS)
end


#Load parameters required for determining initial conditions
ε = 5.3
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
Ace_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_N_0/Yre_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
Are_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_N_0/Yce_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

#In the Southern Region
Ace_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_S_0/Yre_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
Are_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_S_0/Yce_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

# Save de initial state variable in an array
u0= [Are_N_0, Ace_N_0, Are_S_0,Ace_S_0,S_0]

# Save parameters in an array
dt = 5    # 1 Quarterly
D = 200.0 # Simulate for 30 years
N_t = Int(D/dt) # Corresponding no of time steps
index_vector = collect(0:dt:D)

ce_tax_N_int = zeros(length(index_vector))
ce_tax_S_int = zeros(length(index_vector))

# Definimos el modelo de Turing
@model max_molina_ediam() =
begin

    # Define array for the non-state variables
    global N_fossil_energy = []
    global S_fossil_energy = []
    global r_Cost_S_Damage = []
    global r_Consumption_N = []
    global r_Consumption_S = []
    global r_Delta_Temp = []
    global tiempo_line = []
    # Priors
    σ₂ ~ Truncated(Normal(0, 1), 0, Inf)
    ce_tax_S ~ Truncated(Normal(0.2,√σ₂ ),0.1, 0.5)
    ce_tax_N ~ Truncated(Normal(0.6,√σ₂ ),0.3, 0.7)
    Tec_subsidy_N ~ Truncated(Normal(0.2,√σ₂ ),0.1, 0.9)
    RD_subsidy_N ~ Truncated(Normal(0.8,√σ₂ ),0.1, 1.9)
    Tec_subsidy_S ~ Truncated(Normal(0.2,√σ₂ ),0.1, 0.9)
    RD_subsidy_S ~ Truncated(Normal(0.8,√σ₂ ),0.1, 1.9)


    p = [ε ,α ,size_factor,γ_re ,k_re,γ_ce ,k_ce,η_re,η_ce,ν_re ,ν_ce,qsi ,δ_S ,Δ_T_Disaster,β_T ,CO2_base ,CO2_Disaster,labor_growth_N ,labor_growth_S,ρ,λ,σ,
        index_vector,ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S]

    # Solve the ODE system
    prob1 = ODEProblem(molina_ediam,u0,(0.0,D),p)
    global sol = solve(prob1,Euler(),dt=dt)

    # Use simulation output to estimate value of objective function
    utility_consumer_N = 1 .+ ((r_Cost_S_Damage .* r_Consumption_N).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))
    utility_consumer_S = 1 .+ ((r_Cost_S_Damage .* r_Consumption_S).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))

    TOTAL_UTILITY = sum(utility_consumer_N) +sum(utility_consumer_S)
    # Likelihood
    OBJ ~ Truncated(Normal(TOTAL_UTILITY, 10),40,Inf)
    println(OBJ)

end;

# I use optim for the maximum Likelihood estimation
model_optim = max_molina_ediam()
obj_max = 0
mle_estimate = 0
while obj_max < 100
    #mle_estimate = optimize(model_optim,MAP(), NelderMead())
    mle_estimate = optimize(model_optim,MLE(), NelderMead(), Optim.Options(iterations=10_000, allow_f_increases=true))
    obj_max = mle_estimate.values[:OBJ]
end

mc_optim = sample(model_optim, SMC(), MCMCThreads(), 1000, 4,init_theta = mle_estimate.values.array)
plot(mc_optim)

# Greficamos el incremento de la temperatura
plot([2012+i for i in tiempo_line],r_Delta_Temp,
 title = "Incremento de la temperatura",
  lw = 3)


global tmt_agrega = []

for i in 1:30
    println(i)
    mc_optim = sample(model_optim, SMC(), MCMCThreads(), 1000, 4,init_theta = mle_estimate.values.array)
    push!(tmt_agrega,r_Delta_Temp)
end

med_b = mean(hcat(tmt_agrega))
min_b = minimum(hcat(tmt_agrega))
max_b = maximum(hcat(tmt_agrega))
upper_bound = mean(hcat(tmt_agrega)) .+ 1.97 .* std(hcat(tmt_agrega))
lower_bound = mean(hcat(tmt_agrega)) .- 1.97 .* std(hcat(tmt_agrega))

scatter(([2012+i for i in tiempo_line], med_b),
		yerror=(med_b-min_b,
				max_b-med_b))

#=
REFERENCIAS

https://discourse.julialang.org/t/solving-an-optimal-control-problem-with-jump/36570
https://docs.juliahub.com/DiffEqFlux/BdO4p/1.10.3/examples/LV-ODE/
=#
