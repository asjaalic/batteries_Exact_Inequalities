# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       #, state_variables::states When we have 2 hydropower plants- 2 turbines

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc, NSteps) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         

    k = max_SOC/(2*Nfull)
    Small = 0.64

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.01)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps] <= 1, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps] <= 1, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps] <= Small, base_name = "Degradation")

    @variable(M, min_SOH <= soh_final[iStage=1:NStages] <= max_SOH, base_name = "Final_Capacity")        #energy_Capacity     [iStage=1:NStages]
    @variable(M, min_SOH <= soh_new[iStage=1:NStages] <= max_SOH, base_name = "Initial_Capacity")     #energy_Capacity      [iStage=1:NStages]

    #VARIABLES FOR ENVELOPES

    @variable(M, x[iStep=1:NSteps+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps+1], Bin, base_name = "Binary_3")
    @variable(M, u[iStep=1:NSteps+1], Bin, base_name = "Binary_4")

    @variable(M, 0<= w_xx[iStep=1:NSteps+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps+1] <= 1, base_name = "yz")

    @variable(M, 0 <= w_uu[iStep=1:NSteps+1] <=1, base_name = "uu")
    @variable(M, 0 <= w_xu[iStep=1:NSteps+1] <=1, base_name = "xu")
    @variable(M, 0 <= w_yu[iStep=1:NSteps+1] <=1, base_name = "yu")
    @variable(M, 0 <= w_zu[iStep=1:NSteps+1] <=1, base_name = "zu")

    # DEFINE OJECTIVE function - length(Battery_price) = NStages+1=21

    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=1:NSteps) -
      sum(Battery_price[iStage]*(soh_new[iStage]-soh_final[iStage-1]) for iStage=2:NStages) - 
      Battery_price[1]*(soh_new[1]-min_SOH) + 
      Battery_price[NStages+1]*(soh_final[NStages]-min_SOH) 
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M,energy[iStep=1:NSteps], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*max_P*NHoursStep == soc[iStep+1] )

    #@constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
    @constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]) == soc[iStep])    
    
   #@constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
    @constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_zu[iStep]+64*w_uu[iStep])*((max_SOC-min_SOC)/disc)^2)

    # INEQUALITIES CONSTRAINTS
    @constraint(M, xx_1[iStep=1:NSteps+1], w_xx[iStep] <= x[iStep])
    @constraint(M, xx_2[iStep=1:NSteps+1], w_xx[iStep] >= 2*x[iStep]-1)

    @constraint(M, xy_1[iStep=1:NSteps+1], w_xy[iStep] <= x[iStep])
    @constraint(M, xy_2[iStep=1:NSteps+1], w_xy[iStep] <= y[iStep])
    @constraint(M, xy_3[iStep=1:NSteps+1], w_xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M, xz_1[iStep=1:NSteps+1], w_xz[iStep] <= x[iStep])
    @constraint(M, xz_2[iStep=1:NSteps+1], w_xz[iStep] <= z[iStep])
    @constraint(M, xz_3[iStep=1:NSteps+1], w_xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M, yy_1[iStep=1:NSteps+1], w_yy[iStep] <= y[iStep])
    @constraint(M, yy_2[iStep=1:NSteps+1], w_yy[iStep] >= 2*y[iStep]-1)

    @constraint(M, zz_1[iStep=1:NSteps+1], w_zz[iStep] <= z[iStep])
    @constraint(M, zz_2[iStep=1:NSteps+1], w_zz[iStep] >= 2*z[iStep]-1)

    @constraint(M, zy_1[iStep=1:NSteps+1], w_zy[iStep] <= z[iStep])
    @constraint(M, zy_2[iStep=1:NSteps+1], w_zy[iStep] <= y[iStep])
    @constraint(M, zy_3[iStep=1:NSteps+1], w_zy[iStep] >= z[iStep]+y[iStep]-1)

    @constraint(M, uu_1[iStep=1:NSteps+1], w_uu[iStep] <= u[iStep])
    @constraint(M, uu_2[iStep=1:NSteps+1], w_uu[iStep] >= 2*u[iStep]-1)

    @constraint(M, xu_1[iStep=1:NSteps+1], w_xu[iStep] <= x[iStep])
    @constraint(M, xu_2[iStep=1:NSteps+1], w_xu[iStep] <= u[iStep])
    @constraint(M, xu_3[iStep=1:NSteps+1], w_xu[iStep] >= x[iStep]+u[iStep]-1)

    @constraint(M, yu_1[iStep=1:NSteps+1], w_yu[iStep] <= y[iStep])
    @constraint(M, yu_2[iStep=1:NSteps+1], w_yu[iStep] <= u[iStep])
    @constraint(M, yu_3[iStep=1:NSteps+1], w_yu[iStep] >= y[iStep]+u[iStep]-1)

    @constraint(M, zu_1[iStep=1:NSteps+1], w_zu[iStep] <= z[iStep])
    @constraint(M, zu_2[iStep=1:NSteps+1], w_zu[iStep] <= u[iStep])
    @constraint(M, zu_3[iStep=1:NSteps+1], w_zu[iStep] >= z[iStep]+u[iStep]-1)


    # CONSTRAINTS ON DEGRADATION

    @constraint(M, deg_1[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING

    @constraint(M,soh[iStage=1:(NStages-1)], soh_new[iStage+1] >= soh_final[iStage])      #soh[iStage=1:(NStages-1)]

    @constraint(M,final_soh[iStage=1:NStages], soh_final[iStage] == soh_new[iStage] - sum(deg[iStep] for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1]))*k )     #deg2


    return BuildStageProblem(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        deg,
        x,
        y,
        z,
        u,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        w_uu,
        w_xu,
        w_yu,
        w_zu,
        soh_final,
        soh_new,
      )
end

