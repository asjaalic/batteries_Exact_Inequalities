# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       #, state_variables::states When we have 2 hydropower plants- 2 turbines

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NSteps, NHoursStep, NHoursStage, conv, disc) = InputParameters;  
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         

    k = NHoursStep/(2*Nfull)
    #k= NHoursStep*max_P/(2*Nfull*max_SOC)
    #deg_max= (1/Eff_discharge)^0.5        # 1.05
    #deg_avg = deg_max/2                   # 0.21



    M = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(M, "MIPGap", MIPGap)
    #set_optimizer_attribute(M,"MIPFocus", MIPFocus)
    #set_optimizer_attribute(M,"Method", Method)
    #set_optimizer_attribute(M,"Cuts", Cuts)
    #set_optimizer_attribute(M,"Heuristics", Heuristics)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC <= soc_quad[iStep=1:NSteps+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps] <= max_P, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps] <= Big, base_name = "Degradation")

    @variable(M, min_SOH <= soh_final[iStage=1:NStages] <= max_SOH, base_name = "Final_Capacity")        #energy_Capacity
    @variable(M, min_SOH <= soh_new[iStage=1:NStages] <= max_SOH, base_name = "Initial_Capacity")     #energy_Capacity

    #VARIABLES FOR ENVELOPES

    @variable(M, x[iStep=1:NSteps+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps+1], Bin, base_name = "Binary_3")
   # @variable(M, w[iStep=1:NSteps+1], Bin, base_name = "Binary_4")

    @variable(M, 0<= w_xx[iStep=1:NSteps+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps+1] <= 1, base_name = "yz")

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
    @constraint(M, soc[1]== soc[end])

    @constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
 
    @constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
    
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

    # CONSTRAINTS ON DEGRADATION

    @constraint(M, deg_1[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING

    @constraint(M,soh[iStage=1:(NStages-1)], soh_new[iStage+1] >= soh_final[iStage])

    @constraint(M,final_soh[iStage=1:NStages], soh_final[iStage] == soh_new[iStage]- sum(deg[iStep] for iStep=((iStage-1)*NHoursStage+1):(NHoursStage*iStage))*k )     #deg2

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
        #w,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        #SOC_aux,
        #P_aux,
        #d,
        #u,
        #d_1,
        #d_2,
        #deg_1,
        #deg_2,
        soh_final,
        soh_new,
      )
end



  # LINEAR LINEARIZATION
  #=
#CONSTRAINT FOR CHARGING - DEGRADATION
    @constraint(M, deg_pos_1[iStep=1:NSteps], deg2[iStep] >= Eff_charge*0*2+0*Eff_charge*(auxiliary[iStep]-2)+2*Eff_charge*(charge[iStep]-0))
    @constraint(M, deg_pos_2[iStep=1:NSteps], deg2[iStep] >= Eff_charge*0*1.5+0*Eff_charge*(auxiliary[iStep]-1.5)+1.5*Eff_charge*(charge[iStep]-0))
    @constraint(M, deg_pos_3[iStep=1:NSteps], deg2[iStep] >= Eff_charge*0*1+0*Eff_charge*(auxiliary[iStep]-1)+1*Eff_charge*(charge[iStep]-0))
    @constraint(M, deg_pos_4[iStep=1:NSteps], deg2[iStep] >= Eff_charge*0*0.5+0*Eff_charge*(auxiliary[iStep]-0.5)+0.5*Eff_charge*(charge[iStep]-0))

    @constraint(M, deg_pos_5[iStep=1:NSteps], deg2[iStep] >= Eff_charge*2*1.82+2*Eff_charge*(auxiliary[iStep]-1.82)+1.82*Eff_charge*(charge[iStep]-2))
    @constraint(M, deg_pos_6[iStep=1:NSteps], deg2[iStep] >= Eff_charge*2*1.32+2*Eff_charge*(auxiliary[iStep]-1.32)+1.32*Eff_charge*(charge[iStep]-2))
    @constraint(M, deg_pos_7[iStep=1:NSteps], deg2[iStep] >= Eff_charge*2*0.82+2*Eff_charge*(auxiliary[iStep]-0.82)+0.82*Eff_charge*(charge[iStep]-2))
    @constraint(M, deg_pos_8[iStep=1:NSteps], deg2[iStep] >= Eff_charge*2*0.32+2*Eff_charge*(auxiliary[iStep]-0.32)+0.32*Eff_charge*(charge[iStep]-2))

    @constraint(M, deg_pos_9[iStep=1:NSteps], deg2[iStep] >= Eff_charge*4*1.64+4*Eff_charge*(auxiliary[iStep]-1.64)+1.64*Eff_charge*(charge[iStep]-4))
    @constraint(M, deg_pos_10[iStep=1:NSteps], deg2[iStep] >= Eff_charge*4*1.14+4*Eff_charge*(auxiliary[iStep]-1.14)+1.14*Eff_charge*(charge[iStep]-4))
    @constraint(M, deg_pos_11[iStep=1:NSteps], deg2[iStep] >= Eff_charge*4*0.64+4*Eff_charge*(auxiliary[iStep]-0.64)+0.64*Eff_charge*(charge[iStep]-4))

    @constraint(M, deg_pos_12[iStep=1:NSteps], deg2[iStep] >= Eff_charge*6*1.46+6*Eff_charge*(auxiliary[iStep]-1.46)+1.46*Eff_charge*(charge[iStep]-6))
    @constraint(M, deg_pos_13[iStep=1:NSteps], deg2[iStep] >= Eff_charge*6*0.96+6*Eff_charge*(auxiliary[iStep]-0.96)+0.96*Eff_charge*(charge[iStep]-6))

    @constraint(M, deg_pos_14[iStep=1:NSteps], deg2[iStep] >= Eff_charge*8*1.28+8*Eff_charge*(auxiliary[iStep]-1.28)+1.28*Eff_charge*(charge[iStep]-8))
    @constraint(M, deg_pos_15[iStep=1:NSteps], deg2[iStep] >= Eff_charge*8*0.78+8*Eff_charge*(auxiliary[iStep]-0.78)+0.78*Eff_charge*(charge[iStep]-8))

    @constraint(M, deg_pos_16[iStep=1:NSteps], deg2[iStep] >= Eff_charge*10*1.1+10*Eff_charge*(auxiliary[iStep]-1.1)+1.1*Eff_charge*(charge[iStep]-10))
    

    #CONSTRAINTS FOR DISCHARGING - DEGRADATION
    @constraint(M, deg_neg_1[iStep=1:NSteps], deg1[iStep] >= (0*2+0*(auxiliary[iStep]-2)+2*(discharge[iStep]-0))/Eff_discharge)
    @constraint(M, deg_neg_2[iStep=1:NSteps], deg1[iStep] >= (0*1.5+0*(auxiliary[iStep]-1.5)+1.5*(discharge[iStep]-0))/Eff_discharge)
    @constraint(M, deg_neg_3[iStep=1:NSteps], deg1[iStep] >= (0*1+0*(auxiliary[iStep]-1)+1*(discharge[iStep]-0))/Eff_discharge)
    @constraint(M, deg_neg_4[iStep=1:NSteps], deg1[iStep] >= (0*0.5+0*(auxiliary[iStep]-0.5)+0.5*(discharge[iStep]-0))/Eff_discharge)

    @constraint(M, deg_neg_5[iStep=1:NSteps], deg1[iStep] >= (2*1.72+2*(auxiliary[iStep]-1.72)+1.72*(discharge[iStep]-2))/Eff_discharge)
    @constraint(M, deg_neg_6[iStep=1:NSteps], deg1[iStep] >= (2*1.22+2*(auxiliary[iStep]-1.22)+1.22*(discharge[iStep]-2))/Eff_discharge)
    @constraint(M, deg_neg_7[iStep=1:NSteps], deg1[iStep] >= (2*0.72+2*(auxiliary[iStep]-0.72)+0.72*(discharge[iStep]-2))/Eff_discharge)
    @constraint(M, deg_neg_8[iStep=1:NSteps], deg1[iStep] >= (2*0.22+2*(auxiliary[iStep]-0.22)+0.22*(discharge[iStep]-2))/Eff_discharge)

    @constraint(M, deg_neg_9[iStep=1:NSteps], deg1[iStep] >= (4*1.44+4*(auxiliary[iStep]-1.44)+1.44*(discharge[iStep]-4))/Eff_discharge)
    @constraint(M, deg_neg_10[iStep=1:NSteps], deg1[iStep] >= (4*0.94+4*(auxiliary[iStep]-0.94)+0.94*(discharge[iStep]-4))/Eff_discharge)
    @constraint(M, deg_neg_11[iStep=1:NSteps], deg1[iStep] >= (4*0.44+4*(auxiliary[iStep]-0.44)+0.44*(discharge[iStep]-4))/Eff_discharge)

    @constraint(M, deg_neg_12[iStep=1:NSteps], deg1[iStep] >= (6*1.17+6*(auxiliary[iStep]-1.17)+1.17*(discharge[iStep]-6))/Eff_discharge)
    @constraint(M, deg_neg_13[iStep=1:NSteps], deg1[iStep] >= (6*0.67+6*(auxiliary[iStep]-0.67)+0.67*(discharge[iStep]-6))/Eff_discharge)

    @constraint(M, deg_neg_14[iStep=1:NSteps], deg1[iStep] >= (8*0.89+8*(auxiliary[iStep]-0.89)+0.89*(discharge[iStep]-8))/Eff_discharge)

  =#