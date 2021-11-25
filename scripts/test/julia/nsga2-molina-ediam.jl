using ForwardDiff
using Printf
using Random
using LinearAlgebra
using Distributions
using DifferentialEquations
using Plots,StatsPlots
#=
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Algoritmo genético
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
=#

function torneo_k_enario(y,k)
    p = randperm(length(y))
    return p[argmin(y[p[1:k]])]
end

function muestra_uniforme(m,a,b)
    return [a+rand(length(a)).*(b-a) for i in 1:m]
end

function cruza_SBX(p1,p2)
    u = rand()
    nc = 1

    if u <= 0.5
        β = (2*u)^(1/(nc+1))
    else
        β =(1/(2*(1-u)))^(1/(nc+1))
    end

    return [abs.(0.5 .* ((p1+p2) - β .*abs.(p2-p1))),abs.(0.5 .* ((p1+p2) + β .*abs.(p2-p1)))]

end

function mutacion_polinomial(hijo,σ)
    return hijo + randn(length(hijo))*σ
end

function algoritmo_genetico(f,poblacion,seleccion,k,it_max)

    for it in 1:it_max
        @printf("Iteración %d \n", it)
        y = f.(poblacion);
        padres = [[torneo_k_enario(y,k),torneo_k_enario(y,k)] for i in y];
        hijos =  [cruza_SBX(poblacion[p[1]],poblacion[p[2]]) for p in padres];
        poblacion = [mutacion_polinomial(hijo,0.1) for hijo in hijos];
    end

    return poblacion[argmin(f.(poblacion))]
end

dominates(y, y′) = all(y .≤ y′) && any(y .< y′)

function non_dominated_sorting(y)

    Sp_dict = Dict{Int64,Set{Int64}}(i => Set() for i in 1:length(y))
    np_dict = Dict{Int64,Int64}(i => 0 for i in 1:length(y))
    𝐹 = Dict{Int64,Set{Int64}}(1 => Set())
    rank = Dict{Int64,Int64}()

    for (i,p) ∈ enumerate(y)
        for (j,q) ∈ enumerate(y)
            if i !=j
                # Si p domina a q, agregamos q a las soluciones
                # dominadas por p
                if dominates(p,q)
                    push!(Sp_dict[i],j)
                # En caso que q domine a p, incrementamos el contador de soluciones
                # que dominan a p
                elseif dominates(q,p)
                    np_dict[i]+=1
                end
            end
        end
    if np_dict[i]==0
        rank[i]= 1
        push!(𝐹[1],i)
    end
    end

    # Inicializamos el contador de frentes
    i = 1
    while !isempty(𝐹[i])
        𝑄 = Set{Int64}()

        for p ∈ 𝐹[i]
            for q ∈ Sp_dict[p]
                np_dict[q] -=1
                if np_dict[q]==0
                    rank[q] = i + 1
                    push!(𝑄,q)
                end
            end
        end

        i +=1
        𝐹[i] = 𝑄
    end

    return delete!(𝐹,i)
end

function crowding_distance(F,generacion,k)
    valores = [generacion[i] for i in keys(F.dict)]
    max_fi = Dict{Int64,Any}()
    min_fi = Dict{Int64,Any}()

    for i in 1:k
        min_fi[i] = minimum([valor[i] for valor in valores])
        max_fi[i] = maximum([valor[i] for valor in valores])
    end

    # Inicializamos CD de la j-ésima solución en F
    P_CD = Vector{Float64}(undef, length(valores))

    for i in 1:k
        P_sort = sort([(v,i) for (i,v) in enumerate([valor[i] for valor in valores])])
        P_CD[P_sort[1][2]] = Inf16

        P_CD[P_sort[length(valores)][2]] = Inf16

        for j in 2:(length(valores)-1)
            P_CD[P_sort[j][2]]+= ( (P_sort[j+1][1] - P_sort[j-1][1])/(max_fi[i]-min_fi[i]) )
        end
    end

    # [(v,i) for (i,v) in enumerate([valor[i] for valor in P_CD])]
    ordena = [i for (v,i) in sort([(v,i) for (i,v) in enumerate(P_CD)],rev=true)]
    return valores[ordena]

end

function rm_values(poblacion)

    if maximum(poblacion) > 1
        return false
    else
        return true
    end

end

function entre_cero_uno(x)
    if x < 1.0
        nx = x
    else
        sobrante = x-1
        nx  = 1 - sobrante
    end
    return abs(nx)
end

function nsga_II(f,m,x_dim,ktorneo,generaciones,kobjetivos,m_mantiene)
    # Generamos una población aleatoria de tamaño M
    poblacion = muestra_uniforme(m,zeros(x_dim),ones(x_dim));
    # Asignamos un ranking de acuerdo a la dominancia de Pareto
    y = f.(poblacion)
    #non_dominated_sorting(y)
    # Generamos la generación de hijos
    # Torneo binario para elegir a los padres
    padres = [[torneo_k_enario([sum(i) for i in y],ktorneo),torneo_k_enario([sum(i) for i in y],ktorneo)] for i in y];
    # Recombinación y mutación
    hijos =  [i for p in padres for i in cruza_SBX(poblacion[p[1]],poblacion[p[2]])]
    hijos = [entre_cero_uno.(mutacion_polinomial(hijo,0.15)) for hijo in hijos];
    generacion_actual = [i for i in vcat(poblacion,hijos) if rm_values(i)]

    for i ∈ 1:generaciones
        println("Generación $i")
        # Con la población de padres e hijos
        # Asignamos un ranking de acuerdo a la dominancia de Pareto
        # Evaluamos los valores objetivo
        y = f.(generacion_actual);
        #println("Asignamos un ranking de acuerdo a la dominancia de Pareto")
        fronteras = non_dominated_sorting(y)
        # Hacemos un loop y agregamos la solucion a la siguiente generación comenzando
        # del primer frente hasta que tengamos m individuos
        nueva_generacion = []

        i = 1
        while length(nueva_generacion)!=m_mantiene
            #println("Evaluamos en la frontera $i, faltan ", (m_mantiene - length(nueva_generacion)))
            #println( (m - length(nueva_generacion)))
            if length(fronteras[i]) <= (m_mantiene - length(nueva_generacion))
                nueva_generacion = vcat(nueva_generacion,[generacion_actual[i] for i in keys(fronteras[i].dict)]);
                i+=1
            else
                cd = crowding_distance(fronteras[i],generacion_actual,kobjetivos);
                nueva_generacion = vcat(nueva_generacion,cd[1:(m_mantiene - length(nueva_generacion))]);
            end
            #println("Tenemos que llenar más el arreglo?", length(nueva_generacion)!=m_mantiene)
        end

        # Creamos la siguiente generacion con torneo binario, recombinación y mutación
        padres = [[torneo_k_enario(nueva_generacion,ktorneo),torneo_k_enario(nueva_generacion,ktorneo)] for i in nueva_generacion];
        generacion_actual =  [i for p in padres for i in cruza_SBX(nueva_generacion[p[1]],nueva_generacion[p[2]])];
        generacion_actual = [entre_cero_uno.(mutacion_polinomial(hijo,0.15)) for hijo in generacion_actual];
        generacion_actual = [i for i in generacion_actual if rm_values(i)]


    end

    return  generacion_actual
end


function ode_ediam(du,u,p,t)
    ## Parámetros iniciales
    ε,α,size_factor,γ_re,k_re,γ_ce,k_ce,η_re,η_ce,ν_re,ν_ce,qsi,δ_S,Δ_T_Disaster,β_T,CO2_base,CO2_Disaster,labor_growth_N,labor_growth_S,ρ,λ,σ,ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S = p

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

function solve_ediam(X)
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
    Δ_T_Disaster= 6
    labor_growth_N = 0.000
    labor_growth_S = 0.000
    ρ = 0.015
    λ = 0.1443
    σ = 2


    ### Parametros. Global
    β_T,CO2_base,CO2_Disaster,qsi,δ_S ,S_0,Yre_N_0,Yce_N_0 ,Yre_S_0,Yce_S_0  = [3.189799 293.071388 3076.849744 0.010023 0.000630 2695.295397 45.55074 193.2 27.82166 257.5463]

    #Initial Productivity conditions are determined by the initial levels of production of energy
    #In the Northern Region
    global Ace_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_N_0/Yre_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
    global Are_N_0 = ((Yce_N_0^((ε-1)/ε)+Yre_N_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_N_0/Yce_N_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

    #In the Southern Region
    Ace_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yce_S_0/Yre_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))
    Are_S_0 = (1/size_factor)*((Yce_S_0^((ε-1)/ε)+Yre_S_0^((ε-1)/ε))^(ε/(ε-1)))*(1+(Yre_S_0/Yce_S_0)^((1-ε)/ε))^(1/((1-α)*(1-ε)))

    # Save de initial state variable in an array
    u0= [Are_N_0, Ace_N_0, Are_S_0,Ace_S_0,S_0]

    # Save parameters in an array
    dt = 5    # 1 Quarterly
    D = 120.0 # Simulate for 30 years
    N_t = Int(D/dt) # Corresponding no of time steps


    # Define array for the non-state variables
    global N_fossil_energy = []
    global S_fossil_energy = []
    global r_Cost_S_Damage = []
    global r_Consumption_N = []
    global r_Consumption_S = []
    global r_Delta_Temp = []
    global tiempo_line = []

    ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S = X

    p = [ε ,α ,size_factor,γ_re ,k_re,γ_ce ,k_ce,η_re,η_ce,ν_re ,ν_ce,qsi ,δ_S ,Δ_T_Disaster,β_T ,CO2_base ,CO2_Disaster,labor_growth_N ,labor_growth_S,ρ,λ,σ,
        ce_tax_N,ce_tax_S,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S]

    # Solve the ODE system
    prob1 = ODEProblem(ode_ediam,u0,(0.0,D),p)
    global sol = solve(prob1,Euler(),dt=dt)

    # Use simulation output to estimate value of objective function
    utility_consumer_N = 1 .+ ((r_Cost_S_Damage .* r_Consumption_N).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))
    utility_consumer_S = 1 .+ ((r_Cost_S_Damage .* r_Consumption_S).^(1-σ) / (1- σ)) .* (1 ./ ((1+ ρ).^(sol.t)))

    TOTAL_UTILITY = sum(utility_consumer_N) +sum(utility_consumer_S)
    MXTMP = maximum(r_Delta_Temp)

    return TOTAL_UTILITY,-MXTMP
end

sol_nsgaII = nsga_II(solve_ediam,8000,6,2,200,2,200)
A = solve_ediam.(sol_nsgaII);
nds = non_dominated_sorting(A)
A_sub = [A[i] for i in nds[1]]


pareto = scatter([abs(a[2]) for a in A_sub],[abs(a[1]) for a in A_sub], label ="Multi-objetivo",legend=:bottomright, title ="Frontera de Pareto")
scatter!(pareto,[5.321413853203894],[49.83994188094781], label="Mono-objetivo")
