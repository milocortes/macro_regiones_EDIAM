library(deSolve)
library(ggplot2)

Run.ID <- 0
EndTime <- 30
TimeStep <- 1
#Determine policy parameters
#lists of policies
#P0.FWA
#P1.Nordhaus
#P2.Nordhauds+TechnologyPolicy
#P2.1."Nordhaus+TechnologyPolicy.Both"
#P3.Nordhaus+TraditionalGreenClimateFund
#P3.1."Nordhaus+TraditionalGreenClimateFund+R&DS"
#P3.2."Nordhaus+CoR&DGreenClimateFund"
#P3.3."Nordhaus+CoR&DGreenClimateFund+TecS"
#P4.Nordhaus+R&DGreenClimateFund
policy.vector<-c(
  #Start time of policies
  policy.start.time = 0.033,
  #Duration of policies
  #policy.duration=signif(as.numeric(policies[2]), digits = 2),
  #carbon tax
  tax.rate.N=0,
  tax.rate.S=0,
  #Technology push in North
  epsi.re.subsidy.N = 0,

  s.re.subsidy.N = 0,
  #Traditional Green Climate Fund
  epsi.re.subsidy.S = 0,

  # R&D Green Climate Fund
  s.re.subsidy.S = 0,
  s.re.GFsubsidy.N = 0,
  epsi.re.GFsubsidy.N=0)
#Load parameters required for determining initial conditions
alfa <- 0.33
epsilon <- 3.5
size.factor<- 1
Yre.N.0<-25.1
Yce.N.0<-144.9
Yre.S.0<-9.0
Yce.S.0<-105.3
S.0<-915.970085


#Initial Productivity conditions are determined by the initial levels of production of energy
#In the Northern Region
Ace.N.0<-((Yce.N.0^((epsilon-1)/epsilon)+Yre.N.0^((epsilon-1)/epsilon))^(epsilon/(epsilon-1)))*(1+(Yce.N.0/Yre.N.0)^((1-epsilon)/epsilon))^(1/((1-alfa)*(1-epsilon)))
Are.N.0<-((Yce.N.0^((epsilon-1)/epsilon)+Yre.N.0^((epsilon-1)/epsilon))^(epsilon/(epsilon-1)))*(1+(Yre.N.0/Yce.N.0)^((1-epsilon)/epsilon))^(1/((1-alfa)*(1-epsilon)))

#In the Southern Region
Ace.S.0<-(1/size.factor)*((Yce.S.0^((epsilon-1)/epsilon)+Yre.S.0^((epsilon-1)/epsilon))^(epsilon/(epsilon-1)))*(1+(Yce.S.0/Yre.S.0)^((1-epsilon)/epsilon))^(1/((1-alfa)*(1-epsilon)))
Are.S.0<-(1/size.factor)*((Yce.S.0^((epsilon-1)/epsilon)+Yre.S.0^((epsilon-1)/epsilon))^(epsilon/(epsilon-1)))*(1+(Yre.S.0/Yce.S.0)^((1-epsilon)/epsilon))^(1/((1-alfa)*(1-epsilon)))

InitialConditions <- c(Are.N = Are.N.0,
                       Ace.N = Ace.N.0,
                       Are.S = Are.S.0,
                       Ace.S = Ace.S.0,
                       S = S.0)

#Put all parameters together
Parameters <- c(alfa = alfa,
                epsilon = epsilon,
                size.factor= size.factor,
                Gamma.re = 0.25 ,
                k.re = 0,
                Gamma.ce = 0.25,
                k.ce = 0,
                Eta.re= 0.03,
                Eta.ce= 0.021,
                Nu.re = 0.3,
                Nu.ce= 0.25,
                qsi= 0.010054,
                Delta.S =0.001823 ,
                Delta.Temp.Disaster =7.5 ,
                Beta.Delta.Temp = 4.997053,
                CO2.base =  289.415046,
                CO2.Disaster =1298.216153 ,
                labor.growth.N = 0.008,
                labor.growth.S = 0.015,
                rho = 0.01,
                lambda.S = 0.1443,
                sigma.utility = 2,
                policy.vector)

ModelEngine <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    print(Time)
    #Policy Instruments
    #status of policy objectives
    dS.lag<-ifelse(Time-1<=0,0,lagderiv(Time-1,5))
    R.A.North<-ifelse(Time-1<=0,10,lagvalue(Time-1,2)/lagvalue(Time-1,1))
    R.A.South<-ifelse(Time-1<=0,10,lagvalue(Time-1,4)/lagvalue(Time-1,3))

    #policies in the North
    ce.tax.N<-tax.rate.N
    RD.subsidy.N<-s.re.subsidy.N
    RD.subsidy.GF.N<-s.re.GFsubsidy.N
    Tec.subsidy.N<-epsi.re.subsidy.N
    Tec.subsidy.GF.N<-epsi.re.GFsubsidy.N

    #policies in the South
    ce.tax.S<-tax.rate.S
    RD.subsidy.S<-s.re.subsidy.S+RD.subsidy.GF.N
    Tec.subsidy.S<-epsi.re.subsidy.S+Tec.subsidy.GF.N

    #Economic structure
    #Auxiliaries for both regions
    Phi<-(1-alfa)*(1-epsilon)
    epsi.re<-alfa^2 #this is the cost of production of clean technologies
    epsi.ce<-alfa^2 #this is the cost of production of dirty technologies

    #North Region
    #Auxiliaries in North
    L.N<-exp(labor.growth.N*Time)
    Gamma.re.t.N<-Gamma.re*exp(-k.re*(Are.N/Are.N.0-1)) #gamma displays decreasing returns as in Stiligtz
    Gamma.ce.t.N<-Gamma.ce*exp(-k.ce*(Ace.N/Ace.N.0-1)) #gamma displays decreasing returns as in Stiligtz

    #First we determine the equilibrium levels of relative input prices and relative labor
    RelPrice.N<-((Ace.N/Are.N)^(1-alfa))*(((epsi.re*(1-Tec.subsidy.N))/epsi.ce)^alfa)
    RelLabor.N<-((1+ce.tax.N)^epsilon)*((((1-Tec.subsidy.N)*epsi.re)/epsi.ce)^(alfa*(1-epsilon)))*((Are.N/Ace.N)^(-1*Phi))

    #Second we determine the equilibrium conditions for each sector
    #clean sector
    Labor.re.N<-(RelLabor.N*L.N)/(1+RelLabor.N) #based on the assumption that Labor.re.N+Labor.ce.N=L.N
    Price.re.N<-RelPrice.N/(RelPrice.N^(1-epsilon)+(1)^(1-epsilon))^(1/(1-epsilon)) #based on the assumption that  Price.re.N^(1-epsilon)+Price.ce.N^(1-epsilon)=1
    Agg.demand.re.tech.N<-((((alfa^2)*Price.re.N)/((1-Tec.subsidy.N)*epsi.re))^(1/(1-alfa)))*Labor.re.N*Are.N
    Profits.re.N<-(1+RD.subsidy.N)*Eta.re*epsi.re*((1-alfa)/alfa)*Agg.demand.re.tech.N # Expected profits see annex IV. Equilibrium research profits
    Yre.N<-((((alfa^2)*Price.re.N)/((1-Tec.subsidy.N)*epsi.re))^(alfa/(1-alfa)))*Labor.re.N*Are.N

    #dirty sector
    Labor.ce.N<-L.N/(RelLabor.N+1)
    Price.ce.N<-Price.re.N/RelPrice.N
    Agg.demand.ce.tech.N<-((((alfa^2)*Price.ce.N)/(epsi.ce))^(1/(1-alfa)))*Labor.ce.N*Ace.N
    Profits.ce.N<-Eta.ce*epsi.ce*((1-alfa)/alfa)*Agg.demand.ce.tech.N
    Yce.N<-((((alfa^2)*Price.ce.N)/(epsi.ce))^(alfa/(1-alfa)))*Labor.ce.N*Ace.N

    #Total Production
    Y.N<-((Yre.N)^((epsilon-1)/epsilon)+(Yce.N)^((epsilon-1)/epsilon))^(epsilon/(epsilon-1))

    #Allocation of Scientists
    sre.N<-exp(Profits.re.N)/(exp(Profits.ce.N)+exp(Profits.re.N))
    sce.N<-1-sre.N

    #South Region
    #Auxiliaries in South
    L.S<-(exp(labor.growth.S*Time))*size.factor #the population of the South is 4.6 that of the North,
    Gamma.re.t.S<-Gamma.re
    Gamma.ce.t.S<-Gamma.ce

    #First we determine the equilibrium levels of relative input prices and relative labour
    RelPrice.S<-((Ace.S/Are.S)^(1-alfa))*(((epsi.re*(1-Tec.subsidy.S))/epsi.ce)^alfa)
    RelLabor.S<-((1+ce.tax.S)^epsilon)*((((1-Tec.subsidy.S)*epsi.re)/epsi.ce)^(alfa*(1-epsilon)))*((Are.S/Ace.S)^(-1*Phi))

    #Second we determine the equilibrium conditions for each sector
    #clean sector
    Labor.re.S<-(L.S*RelLabor.S)/(RelLabor.S+1) #based on the assumption that Labor.re.S+Labor.ce.S=L.S
    Price.re.S<-RelPrice.S/(RelPrice.S^(1-epsilon)+(1)^(1-epsilon))^(1/(1-epsilon)) #based on the assumption that  Price.re.S^(1-epsilon)+(Price.ce.S)^(1-epsilon)=1
    Agg.demand.re.tech.S<-((((alfa^2)*Price.re.S)/((1-Tec.subsidy.S)*epsi.re))^(1/(1-alfa)))*Labor.re.S*Are.S
    Profits.re.S<-(1+RD.subsidy.S)*Eta.re*epsi.re*((1-alfa)/alfa)*Agg.demand.re.tech.S # Expected profits see annex IV. Equilibrium research profits
    Yre.S<-((((alfa^2)*Price.re.S)/((1-Tec.subsidy.S)*epsi.re))^(alfa/(1-alfa)))*Labor.re.S*Are.S

    #dirty sector
    Labor.ce.S<-L.S/(RelLabor.S+1)
    Price.ce.S<-Price.re.S/RelPrice.S
    Agg.demand.ce.tech.S<-((((alfa^2)*Price.ce.S)/(epsi.ce))^(1/(1-alfa)))*Labor.ce.S*Ace.S
    Profits.ce.S<-Eta.ce*epsi.ce*((1-alfa)/alfa)*Agg.demand.ce.tech.S # Expected profits see annex IV. Equilibrium research profits
    Yce.S<-((((alfa^2)*Price.ce.S)/(epsi.ce))^(alfa/(1-alfa)))*Labor.ce.S*Ace.S

    #Total Production
    Y.S<-((Yre.S)^((epsilon-1)/epsilon)+(Yce.S)^((epsilon-1)/epsilon))^(epsilon/(epsilon-1))

    #Allocation of Scientists
    sre.S<-exp(Profits.re.S)/(exp(Profits.ce.S)+exp(Profits.re.S))
    sce.S<-1-sre.S


    #Changes in Temperature
    Delta.Temp.Disaster<-Delta.Temp.Disaster #increase in temperature at which there is environmental disaster
    CO2.Concentration<-max(CO2.Disaster-S,CO2.base)
    Delta.Temp<-min(Beta.Delta.Temp*log(CO2.Concentration/CO2.base),Delta.Temp.Disaster)

    #Welfare Calculations
    Consumption.N<-Y.N-epsi.re*Agg.demand.re.tech.N-epsi.ce*Agg.demand.ce.tech.N
    Consumption.S<-(Y.S-epsi.re*Agg.demand.re.tech.S-epsi.ce*Agg.demand.ce.tech.S)*(1/size.factor)
    Cost.S.Damage<-((Delta.Temp.Disaster-Delta.Temp)^lambda.S-lambda.S*Delta.Temp.Disaster^(lambda.S-1)*(Delta.Temp.Disaster-Delta.Temp))/((1-lambda.S)*Delta.Temp.Disaster^lambda.S)


    #Budget restrictions
    Budget.function.N<-ce.tax.N*Price.ce.N*Yce.N-
      Tec.subsidy.N*epsi.re*Agg.demand.re.tech.N - #including costs of technology subsidies
      Tec.subsidy.GF.N*epsi.re*Agg.demand.re.tech.S - #- #green climate fund costs of technology subsidies
      RD.subsidy.N*Eta.re*((epsi.re/alfa)-epsi.re)*Agg.demand.re.tech.N - #costs of research subsidies
      RD.subsidy.GF.N*Eta.re*((epsi.re/alfa)-epsi.re)*Agg.demand.re.tech.S #cost of green climate fund R&D subsidies

    Budget.function.S<-ce.tax.S*Price.ce.S*Yce.S-
      Tec.subsidy.S*epsi.re*Agg.demand.re.tech.S - #including costs of technology subsidies
      RD.subsidy.S*Eta.re*((epsi.re/alfa)-epsi.re)*Agg.demand.re.tech.S  #costs of research subsidies

    #State variables

    #Evolution of Productivity North Region
    dAre.N<-Gamma.re.t.N*Eta.re*sre.N*Are.N
    dAce.N<-Gamma.ce.t.N*Eta.ce*sce.N*Ace.N

    #Evolution of Productivity South Region
    dAre.S<-Gamma.re.t.S*Nu.re*sre.S*(Are.N-Are.S)
    dAce.S<-Gamma.ce.t.S*Nu.ce*sce.S*(Ace.N-Ace.S)

    #Environmental Quality
    dS<-min(1.0,Delta.S*S-qsi*(Yce.N+Yce.S))

    #Define variables to output
    vars.out<-list(c(dAre.N,dAce.N,dAre.S,dAce.S,dS),
                   RelPrice_N = RelPrice.N,
                   RelLabor_N = RelLabor.N,
                   Labor.re_N = Labor.re.N,
                   Price.re_N = Price.re.N,
                   Agg.demand.re.tech_N = Agg.demand.re.tech.N,
                   Profits.re_N = Profits.re.N,
                   Yre_N = Yre.N,
                   sre_N = sre.N,
                   Labor.ce_N = Labor.ce.N,
                   Price.ce_N = Price.ce.N,
                   Agg.demand.ce.tech_N = Agg.demand.ce.tech.N,
                   Profits.ce_N = Profits.ce.N,
                   Yce_N = Yce.N,
                   sce_N = sce.N,
                   Y_N = Y.N,
                   RelPrice_S = RelPrice.S,
                   RelLabor_S = RelLabor.S,
                   Labor.re_S = Labor.re.S,
                   Price.re_S = Price.re.S,
                   Agg.demand.re.tech_S = Agg.demand.re.tech.S,
                   Profits.re_S = Profits.re.S,
                   Yre_S = Yre.S,
                   sre_S = sre.S,
                   Labor.ce_S = Labor.ce.S,
                   Price.ce_S = Price.ce.S,
                   Agg.demand.ce.tech_S = Agg.demand.ce.tech.S,
                   Profits.ce_S = Profits.ce.S,
                   Yce_S = Yce.S,
                   sce_S = sce.S,
                   Y_S = Y.S,
                   Delta.Temp = Delta.Temp,
                   Gamma.re.t_N = Gamma.re.t.N,
                   Gamma.re.t_S = Gamma.re.t.S,
                   Gamma.ce.t_N = Gamma.ce.t.N,
                   Gamma.ce.t_S = Gamma.ce.t.S,
                   L_N = L.N,
                   L_S = L.S,
                   Consumption_N = Consumption.N,
                   Consumption_S = Consumption.S,
                   CO2.Concentration = CO2.Concentration,
                   Cost.S.Damage=Cost.S.Damage,
                   #Policy.Duration = policy.duration*EndTime,
                   Policy.Start.Time = policy.start.time*EndTime,
                   ce.tax_N=ce.tax.N,
                   RD.subsidy_N=RD.subsidy.N,
                   RD.subsidy.GF_N=RD.subsidy.GF.N,
                   Tec.subsidy_N=Tec.subsidy.N,
                   Tec.subsidy.GF_N=Tec.subsidy.GF.N,
                   ce.tax_S=ce.tax.S,
                   RD.subsidy_S=(RD.subsidy.S-RD.subsidy.GF.N),
                   Tec.subsidy_S=(Tec.subsidy.S-Tec.subsidy.GF.N),
                   dS.lag=dS.lag,
                   Budget.function.N = Budget.function.N,
                   Budget.function.S = Budget.function.S)
    return(vars.out)
  })
}


out <- as.data.frame(dede(InitialConditions, seq(0, EndTime, by = TimeStep), ModelEngine, Parameters))
out$Run.ID<-Run.ID
out$time<-out$time+1983

g <- ggplot(out, aes(time))
g <- g + geom_line(aes(y=Yce_S), colour="red")
g <- g + geom_line(aes(y=Yce_N), colour="green")
g

sigma.utility = 2
rho = 0.01


out$Utility.Consumer.N<-ifelse(out$time<(EndTime*0.83),
                               ifelse(out$Cost.S.Damage==0,-100,(1+((out$Cost.S.Damage*out$Consumption_N)^(1-sigma.utility))/(1-sigma.utility))*(1/((1+rho)^out$time))),
                               ifelse(out$Delta.Temp>=2.0,-100,(1+((out$Cost.S.Damage*out$Consumption_N)^(1-sigma.utility))/(1-sigma.utility))*(1/((1+rho)^out$time))))
out$Utility.Consumer.S<-ifelse(out$time<(EndTime*0.83),
                               ifelse(out$Cost.S.Damage==0,-100,(1+((out$Cost.S.Damage*out$Consumption_S)^(1-sigma.utility))/(1-sigma.utility))*(1/((1+rho)^out$time))),
                               ifelse(out$Delta.Temp>=2.0,-100,(1+((out$Cost.S.Damage*out$Consumption_S)^(1-sigma.utility))/(1-sigma.utility))*(1/((1+rho)^out$time))))

#Budget constraint
Budget.N<-sum(out$Budget.function.N)
Budget.S<-sum(out$Budget.function.S)
ifelse(Budget.N<0,-10000,
       ifelse(Budget.S<0,-10000,
              1*(sum(as.numeric(out$Utility.Consumer_N))+sum(as.numeric(out$Utility.Consumer_S)))
       ))
