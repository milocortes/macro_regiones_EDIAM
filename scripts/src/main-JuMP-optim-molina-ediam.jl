include("utils-JuMP.jl")

using DataFrames
using CSV
using Printf
using ProgressBars

# Cargamos los parámetros de las regiones

path = "/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/experimental_design/"

region_parameters_CESM2 = Dict(
                            "america" => CSV.read(path*"cmpi6_CESM2_america_experimental_design.csv",DataFrame),
                            "eurafrica" => CSV.read(path*"cmpi6_CESM2_eurafrica_experimental_design.csv",DataFrame),
                            "asia_oceania" => CSV.read(path*"cmpi6_CESM2_asia_experimental_design.csv",DataFrame)
                        )

region_parameters_GFDL = Dict(
                            "america" => CSV.read(path*"cmpi6_GFDL_america_experimental_design.csv",DataFrame),
                            "eurafrica" => CSV.read(path*"cmpi6_GFDL_eurafrica_experimental_design.csv",DataFrame),
                            "asia_oceania" => CSV.read(path*"cmpi6_GFDL_asia_experimental_design.csv",DataFrame)
                        )

region_parameters = Dict("CESM2" => region_parameters_CESM2, "GFDL" =>region_parameters_GFDL)

policies = ["P"*string(i) for i in 0:7]

gcms = ["CESM2" ,"GFDL"]

resultados = []
for gcm in gcms

    for region in ["america","eurafrica","asia_oceania"]
        printstyled("*******************************************\n"; color = :white)
        printstyled("*******************************************\n"; color = :white)
        printstyled("%%%%%%%%%%%%%  "*region*"  %%%%%%%%%%%%%\n"; color = :yellow)
        printstyled("*******************************************\n"; color = :white)
        printstyled("*******************************************\n"; color = :white)
        println()

        for  parameter_set in ProgressBar(1:400)
            printstyled("++++++++++++++++++++++++++++++++++++++++++\n"; color = :white)
            printstyled("    Parameter set: "*string(parameter_set)*" .Region: "*region*" \n"; color = :yellow)
            printstyled("++++++++++++++++++++++++++++++++++++++++++\n"; color = :white)
            println()

            for policy in policies
                printstyled("-------------------------------------------\n"; color = :white)
                printstyled("    Policy "*policy*" \n"; color = :yellow)
                printstyled("------------------------------------------\n"; color = :white)

                #Load parameters required for determining initial conditions
                ε = region_parameters[gcm][region].epsilon[parameter_set]
                α = region_parameters[gcm][region].alfa[parameter_set]
                size_factor = region_parameters[gcm][region].size_factor[parameter_set]
                γ_re = region_parameters[gcm][region].Gamma_re[parameter_set]
                k_re = region_parameters[gcm][region].k_re[parameter_set]
                γ_ce = region_parameters[gcm][region].Gamma_ce[parameter_set]
                k_ce = region_parameters[gcm][region].k_ce[parameter_set]
                η_re= region_parameters[gcm][region].Eta_re[parameter_set]
                η_ce= region_parameters[gcm][region].Eta_ce[parameter_set]
                ν_re = region_parameters[gcm][region].Nu_re[parameter_set]
                ν_ce= region_parameters[gcm][region].Nu_ce[parameter_set]
                qsi = region_parameters[gcm][region].qsi[parameter_set]
                δ_S = region_parameters[gcm][region].δ_S[parameter_set]
                Δ_T_Disaster = 6
                β_T = region_parameters[gcm][region].β_T[parameter_set]
                CO2_base = region_parameters[gcm][region].CO2_base[parameter_set]
                CO2_Disaster = region_parameters[gcm][region].CO2_Disaster[parameter_set]
                labor_growth_N = region_parameters[gcm][region].labor_growth_N[parameter_set]
                labor_growth_S = region_parameters[gcm][region].labor_growth_S[parameter_set]
                ρ = region_parameters[gcm][region].rho[parameter_set]
                λ = region_parameters[gcm][region].lambda_S[parameter_set]
                σ = region_parameters[gcm][region].sigma_utility[parameter_set]

                ## Y renewable energy, advanced economies
                Yre_N_0 = region_parameters[gcm][region].Yre_N_0[parameter_set]
                ## Y carbon energy, advanced economies
                Yce_N_0 = region_parameters[gcm][region].Yce_N_0[parameter_set]
                ## Y renewable energy, emerging economies
                Yre_S_0 = region_parameters[gcm][region].Yre_S_0[parameter_set]
                ## Y carbon energy, emerging economies
                Yce_S_0 = region_parameters[gcm][region].Yce_S_0[parameter_set]
                ### Environment quality
                S_0 = region_parameters[gcm][region].S_0[parameter_set]

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

                try
                    model = optimize_model(model,policy);

                    solutions = solution_summary(model, verbose=true)

                    ce_tax_S = solutions.primal_solution["ce_tax_S"]
                    ce_tax_N = solutions.primal_solution["ce_tax_N"]
                    Tec_subsidy_N = solutions.primal_solution["Tec_subsidy_N"]
                    RD_subsidy_N = solutions.primal_solution["RD_subsidy_N"]
                    Tec_subsidy_S = solutions.primal_solution["Tec_subsidy_S"]
                    RD_subsidy_S = solutions.primal_solution["RD_subsidy_S"]
                    Tec_subsidy_GF_N = solutions.primal_solution["Tec_subsidy_GF_N"]
                    RD_subsidy_GF_N = solutions.primal_solution["RD_subsidy_GF_N"]

                    optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N)

                    # Verificamos si el incremento de la temperatura sobre pasa los tres grados centígrados
                    println(any(r_Delta_Temp .> 2.0))
                    println(!any(r_Delta_Temp .> 2.0))
                    println(sum(r_Delta_Temp))
                    push!(resultados,!any(r_Delta_Temp .> 2.0))
                catch
                    printstyled("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"; color = :red)
                    printstyled("    ERROR \n"; color = :red)
                    printstyled("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"; color = :red)
                    push!(resultados,"RNF")

                end
            end
        end
    end
end
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

solutions = solution_summary(model, verbose=true)

ce_tax_S = solutions.primal_solution["ce_tax_S"]
ce_tax_N = solutions.primal_solution["ce_tax_N"]
Tec_subsidy_N = solutions.primal_solution["Tec_subsidy_N"]
RD_subsidy_N = solutions.primal_solution["RD_subsidy_N"]
Tec_subsidy_S = solutions.primal_solution["Tec_subsidy_S"]
RD_subsidy_S = solutions.primal_solution["RD_subsidy_S"]
Tec_subsidy_GF_N = solutions.primal_solution["Tec_subsidy_GF_N"]
RD_subsidy_GF_N = solutions.primal_solution["RD_subsidy_GF_N"]

optim_welfare(ce_tax_S,ce_tax_N,Tec_subsidy_N,RD_subsidy_N,Tec_subsidy_S,RD_subsidy_S,Tec_subsidy_GF_N,RD_subsidy_GF_N)
sum(r_Delta_Temp)
# Verificamos si el incremento de la temperatura sobre pasa los tres grados centígrados
any(r_Delta_Temp .> 3.0)
