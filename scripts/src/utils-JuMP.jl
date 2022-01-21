using DifferentialEquations
using Plots,StatsPlots
using JuMP
using Ipopt

function molina_ediam(du,u,p,t)
    ## Parámetros iniciales
    ε,α,size_factor,γ_re,k_re,γ_ce,k_ce,η_re,η_ce,ν_re,ν_ce,qsi,δ_S,Δ_T_Disaster,β_T,CO2_base,CO2_Disaster,labor_growth_N,labor_growth_S,ρ,λ,σ,ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N = p

    #=
    Derivada con respecto al tiempo del vector de estado.
        * u es el vector de estado (arreglo)
    =#
    time =  t

    Are_N,Ace_N,Are_S,Ace_S,S = u

    ### Auxiliares generales

    φ= (1-α)*(1-ε)

    #this is the cost of production of clean technologies
    epsi_re = α^2
    #this is the cost of production of dirty technologies
    epsi_ce = α^2

    ### North Region
    #Auxiliaries in North

    L_N = Float64(ℯ)^(labor_growth_N*time)

    #gamma displays decreasing returns as in Stiligtz
    γ_re_t_N = γ_re*Float64(ℯ)^(-k_re*(Are_N/Are_N_0-1))

    #gamma displays decreasing returns as in Stiligtz
    γ_ce_t_N = γ_ce*Float64(ℯ)^(-k_ce*(Ace_N/Ace_N_0-1))

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

    sre_N = Float64(ℯ)^(Profits_re_N)/(Float64(ℯ)^(Profits_ce_N)+Float64(ℯ)^(Profits_re_N))
    sce_N = 1-sre_N

    #Auxiliaries in South
    #the population of the South is 4.6 that of the North,
    L_S = (Float64(ℯ)^(labor_growth_S*time))*size_factor
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
    sre_S = Float64(ℯ)^(Profits_re_S)/(Float64(ℯ)^(Profits_ce_S)+Float64(ℯ)^(Profits_re_S))
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
    #Tec_subsidy_GF_N = 0
    #RD_subsidy_GF_N = 0

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

function optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N)
    # Define array for the non-state variables
    global N_fossil_energy = []
    global S_fossil_energy = []
    global r_Cost_S_Damage = []
    global r_Consumption_N = []
    global r_Consumption_S = []
    global r_Delta_Temp = []
    global tiempo_line = []


    p = [ε ,α ,size_factor,γ_re ,k_re,γ_ce ,k_ce,η_re,η_ce,ν_re ,ν_ce,qsi ,δ_S ,Δ_T_Disaster,β_T ,CO2_base ,CO2_Disaster,labor_growth_N ,labor_growth_S,ρ,λ,σ,ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N]
    # Solve the ODE system
    prob1 = ODEProblem(molina_ediam,u0,(0.0,D),p)
    global sol = DifferentialEquations.solve(prob1,Euler(),dt=dt)

    if  0 in r_Cost_S_Damage
        total_utility = -1000.0
    else
        # Use simulation output to estimate value of objective function
        utility_consumer_N = 1 .+ ((r_Cost_S_Damage .* r_Consumption_N).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))
        utility_consumer_S = 1 .+ ((r_Cost_S_Damage .* r_Consumption_S).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))

        total_utility = sum(utility_consumer_N) +sum(utility_consumer_S)
    end

    return total_utility
end

function optimize_model(model,policy_regimen)
    if policy_regimen=="P0"
        #= FWA: Future Without Action
            τᴬ,τᴱ,hᴬ,hᴱ,qᴬ,qᴱ,hᴳ,qᴳ = 0
        =#
        @variable(model,ce_tax_N)
        @NLconstraint(model,ce_tax_N == 0)

        @variable(model,ce_tax_S)
        @NLconstraint(model,ce_tax_S == 0)

        @variable(model,Tec_subsidy_N)
        @NLconstraint(model,Tec_subsidy_N == 0)

        @variable(model,Tec_subsidy_S)
        @NLconstraint(model,Tec_subsidy_S == 0)

        @variable(model,RD_subsidy_N)
        @NLconstraint(model,RD_subsidy_N == 0)

        @variable(model,RD_subsidy_S)
        @NLconstraint(model,RD_subsidy_S == 0)

        @variable(model,Tec_subsidy_GF_N)
        @NLconstraint(model,Tec_subsidy_GF_N == 0)

        @variable(model,RD_subsidy_GF_N)
        @NLconstraint(model, RD_subsidy_GF_N == 0)


        register(model, :optim_welfare, 8, optim_welfare, autodiff=true)
        @NLobjective(model, Max, optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N))
        set_optimizer(model, Ipopt.Optimizer)
        optimize!(model)

        return model

    elseif policy_regimen=="P1"
        #= Independent Carbon Tax [Both]
            τᴬ,τᴱ > 0
            hᴬ,hᴱ,qᴬ,qᴱ,hᴳ,qᴳ = 0
        =#
        @variable(model,ce_tax_N)
        @NLconstraint(model, 0.1 <= ce_tax_N)

        @variable(model,ce_tax_S)
        @NLconstraint(model, 0.1 <= ce_tax_S )

        @variable(model,Tec_subsidy_N)
        @NLconstraint(model,Tec_subsidy_N == 0)

        @variable(model,Tec_subsidy_S)
        @NLconstraint(model,Tec_subsidy_S == 0)

        @variable(model,RD_subsidy_N)
        @NLconstraint(model,RD_subsidy_N == 0)

        @variable(model,RD_subsidy_S)
        @NLconstraint(model,RD_subsidy_S == 0)

        @variable(model,Tec_subsidy_GF_N)
        @NLconstraint(model,Tec_subsidy_GF_N == 0)

        @variable(model,RD_subsidy_GF_N)
        @NLconstraint(model, RD_subsidy_GF_N == 0)


        register(model, :optim_welfare, 8, optim_welfare, autodiff=true)
        @NLobjective(model, Max, optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N))
        set_optimizer(model, Ipopt.Optimizer)
        optimize!(model)

        return model

    elseif policy_regimen=="P2"
        #= Independent Carbon Tax [Both] + Independent Technology-R&D[Both]
        τᴬ,τᴱ,hᴬ,hᴱ,qᴬ,qᴱ > 0
        hᴳ,qᴳ = 0
        =#

        @variable(model,ce_tax_N)
        @NLconstraint(model, 0.1 <= ce_tax_N)

        @variable(model,ce_tax_S)
        @NLconstraint(model, 0.1 <= ce_tax_S )

        @variable(model,Tec_subsidy_N)
        @NLconstraint(model, 0.1 <= Tec_subsidy_N )

        @variable(model,Tec_subsidy_S)
        @NLconstraint(model, 0.1 <= Tec_subsidy_S)

        @variable(model,RD_subsidy_N)
        @NLconstraint(model, 0.1 <= RD_subsidy_N )

        @variable(model,RD_subsidy_S)
        @NLconstraint(model, 0.1 <= RD_subsidy_S )

        @variable(model,Tec_subsidy_GF_N)
        @NLconstraint(model, Tec_subsidy_GF_N == 0 )

        @variable(model,RD_subsidy_GF_N)
        @NLconstraint(model, RD_subsidy_GF_N == 0 )


        register(model, :optim_welfare, 8, optim_welfare, autodiff=true)
        @NLobjective(model, Max, optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N))
        set_optimizer(model, Ipopt.Optimizer)
        optimize!(model)

        return model

    elseif policy_regimen=="P3"
        #= Harmonized Carbon Tax [Both] + Cooperative Technology[GCF]-R&D[AR]
            * No R&D subsidies in emerging region
            τᴬ = τᴱ > 0
            hᴬ , qᴬ > 0
            hᴱ = hᴳ > 0
            qᴱ = qᴳ = 0
        =#

        @variable(model,ce_tax_N)
        @NLconstraint(model, 0.1 <= ce_tax_N)

        @variable(model,ce_tax_S)
        @NLconstraint(model, ce_tax_N == ce_tax_S )

        @variable(model,Tec_subsidy_N)
        @NLconstraint(model, 0.1 <= Tec_subsidy_N )

        @variable(model,Tec_subsidy_S)
        @NLconstraint(model, 0.1 <= Tec_subsidy_S)

        @variable(model,RD_subsidy_N)
        @NLconstraint(model, 0.1 <= RD_subsidy_N )

        @variable(model,RD_subsidy_S)
        @NLconstraint(model, 0.1 <= RD_subsidy_S )

        @variable(model,Tec_subsidy_GF_N)
        @NLconstraint(model, Tec_subsidy_GF_N == 0 )

        @variable(model,RD_subsidy_GF_N)
        @NLconstraint(model, RD_subsidy_GF_N == 0 )


        register(model, :optim_welfare, 8, optim_welfare, autodiff=true)
        @NLobjective(model, Max, optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N))
        set_optimizer(model, Ipopt.Optimizer)
        optimize!(model)

        return model

    elseif policy_regimen=="P4"
        #= Harmonized Carbon Tax [Both] + Cooperative Technology[GCF]- Independent R&D[Both]
            * Independent R&D subsidies
            τᴬ = τᴱ > 0
            hᴬ , qᴬ , qᴱ > 0
            hᴱ = hᴳ > 0
            qᴳ = 0
        =#
    elseif policy_regimen=="P5"
        #= Harmonized Carbon Tax [Both] + Cooperative R&D[GCF]- Independent Technology[AR]
            * No technology subsidies in emerging region
            τᴬ = τᴱ > 0
            hᴬ , qᴬ > 0
            hᴱ = hᴳ = 0
            qᴱ = qᴳ > 0
        =#
    elseif policy_regimen=="P6"
        #= Harmonized Carbon Tax [Both] + Cooperative R&D[GCF]- Independent Technology[Both]
            * Independent technology subsidies in emerging region
            τᴬ = τᴱ > 0
            hᴬ , hᴱ , qᴬ > 0
            qᴱ = qᴳ > 0
            hᴳ = 0
        =#
    elseif policy_regimen=="P7"
        #= Harmonized Carbon Tax  +  Cooperative Technology - R&D[GCF]
            * Independent technology subsidies in emerging region
            τᴬ = τᴱ > 0
            hᴬ , qᴬ > 0
            qᴱ = qᴳ > 0
            hᴱ = hᴳ > 0
        =#
        @variable(model,ce_tax_N)
        @NLconstraint(model, 0.1 <= ce_tax_N)

        @variable(model,ce_tax_S)
        @NLconstraint(model, 0.1 <= ce_tax_S )

        @variable(model,Tec_subsidy_N)
        @NLconstraint(model, 0.1 <= Tec_subsidy_N )

        @variable(model,Tec_subsidy_S)
        @NLconstraint(model, 0.1 <= Tec_subsidy_S)

        @variable(model,RD_subsidy_N)
        @NLconstraint(model, 0.1 <= RD_subsidy_N )

        @variable(model,RD_subsidy_S)
        @NLconstraint(model, 0.1 <= RD_subsidy_S )


        register(model, :optim_welfare, 6, optim_welfare, autodiff=true)
        @NLobjective(model, Max, optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S))
        set_optimizer(model, Ipopt.Optimizer)
        optimize!(model)

        return model

    end
end
