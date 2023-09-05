# EXCEL SAVINGS
#using DataFrames
#using XLSX

function data_saving(InputParameters::InputParam,ResultsOpt::Results)

    @unpack (NYears, NMonths, NStages, NSteps, Big, NHoursStep, NHoursStage, disc) = InputParameters;
    
   # @unpack (charge,discharge, soc,soc_quad, soh_final, soh_initial, revenues_per_stage, deg, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy, gain_stage, cost_rev, deg_stage) = ResultsOpt;  
   @unpack (charge,discharge, soc,soc_quad, soh_final, soh_initial, revenues_per_stage, deg, x, y, z, u, w_xx, w_yy, w_zz, w_uu, w_xy, w_xz, w_zy, w_xu, w_yu, w_zu, gain_stage, cost_rev, deg_stage) = ResultsOpt;
   @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "$NStages stages-max_SOH $max_SOH-min_SOC $min_SOC-Ncycles $Nfull-decreasing-15 livelli"
    nameFile="Final results $a" 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    battery_costs= DataFrame()
    
    general[!, "Stage"] = 1:1:NStages
    general[!,"SOH_initial"] = soh_initial
    general[!,"SOH_final"] = soh_final
    general[!,"Degradation"] = deg_stage
    general[!,"Net_Revenues"] = revenues_per_stage
    general[!,"Gain charge/discharge"] = gain_stage
    general[!,"Cost revamping"] = cost_rer

    battery_costs[!,"Costs €/MWh"] = Battery_price[1:NStages+1]

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = 1:1:NSteps                                      #((iStage-1)*NHoursStage+1):(NHoursStage*iStage)
        steps[!, "Energy_prices €/MWh"] = Power_prices[:]
        steps[!, "SOC MWh"] = soc[:]
        steps[!, "Charge MW"] = charge[:]
        steps[!, "Discharge MW"] = discharge[:]
        steps[!, "SOC_quad MW"] = soc_quad[:]
        steps[!, "Deg -"] = deg[:]
        steps[!, "X"] = x[:]
        steps[!, "Y"] = y[:]
        steps[!, "Z"] = z[:]
        steps[!, "U"] = u[:]
        steps[!, "XX"] = w_xx[:]
        steps[!, "YY"] = w_yy[:]
        steps[!, "ZZ"] = w_zz[:]
        steps[!, "UU"] = w_uu[:]
        steps[!, "XY"] = w_xy[:]
        steps[!, "XZ"] = w_xz[:]
        steps[!, "ZY"] = w_zy[:]
        steps[!, "XU"] = w_xu[:]
        steps[!, "YU"] = w_yu[:]
        steps[!, "ZU"] = w_zu[:]

        XLSX.writetable("$iStage stage $a.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end






