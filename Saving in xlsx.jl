# EXCEL SAVINGS
#using DataFrames
#using XLSX

function data_saving(InputParameters::InputParam,ResultsOpt::Results)

    @unpack (NYears, NMonths, NStages, NSteps, Big, NHoursStep, NHoursStage, disc) = InputParameters;
    #@unpack (charge,discharge, soc, soh_final, soh_initial, revenues_per_stage, deg, soc_aux, p_aux, d, deg, d_1,d_2,d_3, deg_1,deg_2,deg_3,u1,u2, gain_stage, cost_rev, deg_stage) = ResultsOpt;        #cum_energy,
    #@unpack (charge,discharge, soc, soh_final, soh_initial, revenues_per_stage, deg, soc_aux, p_aux, d, deg, d_1,d_2, deg_1,deg_2,u, gain_stage, cost_rev, deg_stage) = ResultsOpt;        #cum_energy,

    @unpack (charge,discharge, soc,soc_quad, soh_final, soh_initial, revenues_per_stage, deg, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy, gain_stage, cost_rev, deg_stage) = ResultsOpt;  
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "$NStages stages - versione prof"
    nameFile="Final results decreasing prices $nameF" 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    battery_costs= DataFrame()
    
    general[!,"SOH_initial"] = soh_initial[:]
    general[!,"SOH_final"] = soh_final[:]
    general[!,"Degradation"] = deg_stage[:]
    general[!,"Net_Revenues"] = revenues_per_stage[:]
    general[!,"Gain charge/discharge"] = gain_stage[:]
    general[!,"Cost revamping"] = cost_rev[:]

    battery_costs[!,"Costs €/MWh"] = Battery_price[:]

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = ((iStage-1)*NHoursStage+1):(NHoursStage*iStage)
        steps[!, "Energy_prices €/MWh"] = Power_prices[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!,"SOC MWh"] = soc[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "Charge MW"] = charge[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "Discharge MW"] = discharge[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "SOC_quad MW"] = soc_quad[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "Deg MWh"] = deg[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "X"] = x[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "Y"] = y[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "Z"] = z[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "XX"] = w_xx[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "YY"] = w_yy[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "ZZ"] = w_zz[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "XY"] = w_xy[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "XZ"] = w_xz[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]
        steps[!, "ZY"] = w_zy[((iStage-1)*NHoursStage+1):(NHoursStage*iStage)]

        XLSX.writetable("$iStage stage $nameF.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end






