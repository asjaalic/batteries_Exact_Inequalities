# EXCEL SAVINGS
#using DataFrames
#using XLSX

function data_saving(InputParameters::InputParam,ResultsOpt::Results)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, NSteps) = InputParameters;       #NSteps,NHoursStage
    
   #@unpack (charge,discharge, soc, revenues_per_stage, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy) = ResultsOpt;  
   @unpack (charge,discharge, soc,soc_quad,deg, revenues_per_stage, x, y, z, u, w_xx, w_yy, w_zz, w_uu, w_xy, w_xz, w_zy, w_xu, w_yu, w_zu,soh_initial,soh_final,deg_stage,revenues_per_stage,gain_stage, cost_rev) = ResultsOpt;
   @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "$NSteps steps - $Eff_charge efficiency quad 06.12"
    nameFile="Final results $a" 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    battery_costs= DataFrame()
    
    general[!,"Stage"] = 1:1:NStages
    general[!,"SOH_initial"] = soh_initial[:]
    general[!,"SOH_final"] = soh_final[:]
    general[!,"Degradation"] = deg_stage[:]
    general[!,"Net_Revenues"] = revenues_per_stage[:]
    general[!,"Gain charge/discharge"] = gain_stage[:]
    general[!,"Cost revamping"] = cost_rev[:]

    battery_costs[!,"Costs €/MWh"] = Battery_price[1:NStages+1]

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = (Steps_stages[iStage]+1):(Steps_stages[iStage+1])
        steps[!, "Energy_prices €/MWh"] = Power_prices[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC MWh"] = soc[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Charge MW"] = charge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Discharge MW"] = discharge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC_quad MW"] = soc_quad[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Deg -"] = deg[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "X"] = x[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Y"] = y[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Z"] = z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "U"] = u[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XX"] = w_xx[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "YY"] = w_yy[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "ZZ"] = w_zz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "UU"] = w_uu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XY"] = w_xy[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XZ"] = w_xz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "ZY"] = w_zy[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XU"] = w_xu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "YU"] = w_yu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "ZU"] = w_zu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]

        XLSX.writetable("$iStage stage $a.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end






