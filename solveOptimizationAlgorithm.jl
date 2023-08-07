# SOLVE OPTIMIZATION PROBLEM

function solveOptimizationProblem(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)

    @unpack (NYears, NMonths, NStages, NSteps, Big, NHoursStep, NHoursStage, disc) = InputParameters;
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull ) = Battery;

    println("Solving Optimization Problem")

    k = NHoursStep/(2*Nfull)

    objective = 0
    revenues_per_stage = zeros(NStages)
    gain_stage = zeros(NStages)
    cost_rev = zeros(NStages)
    deg_stage = zeros(NStages)

    charge = zeros(NSteps)
    discharge = zeros(NSteps)
    soc = zeros(NSteps+1)
    deg = zeros(NSteps)

    soc_quad = zeros(NSteps+1)
    x = zeros(NSteps+1)
    y = zeros(NSteps+1)
    z = zeros(NSteps+1)
    #w = zeros(NSteps+1)
    w_xx = zeros(NSteps+1)
    w_yy = zeros(NSteps+1)
    w_zz = zeros(NSteps+1)
    w_xy = zeros(NSteps+1)
    w_xz = zeros(NSteps+1)
    w_zy = zeros(NSteps+1)

    soh_final = zeros(NStages)
    soh_initial = zeros(NStages)
    

    problem = BuildStageProblem(InputParameters, SolverParameters, Battery)

    #unset_time_limit_sec(problem)
    @timeit to "Solve optimization" optimize!(problem.M)

    if termination_status(problem.M) != MOI.OPTIMAL
        println("NOT OPTIMAL: ", termination_status(problem.M))
    else
        println("Optimization finished")
    end

    @timeit to "Collecting results" begin
        objective = JuMP.objective_value(problem.M)
        
        for iStep=1:NSteps
            soc[iStep] = JuMP.value(problem.soc[iStep])
            charge[iStep] = JuMP.value(problem.charge[iStep])
            discharge[iStep] = JuMP.value(problem.discharge[iStep])
            deg[iStep] = JuMP.value(problem.deg[iStep])

            soc_quad[iStep] = JuMP.value(problem.soc_quad[iStep])
            x[iStep] = JuMP.value(problem.x[iStep])
            y[iStep] = JuMP.value(problem.y[iStep])
            z[iStep] = JuMP.value(problem.z[iStep])
            #w = JuMP.value(problem.w[iStep])
            w_xx[iStep] = JuMP.value(problem.w_xx[iStep])
            w_yy[iStep] = JuMP.value(problem.w_yy[iStep])
            w_zz[iStep] = JuMP.value(problem.w_zz[iStep])
            w_xy[iStep] = JuMP.value(problem.w_xy[iStep])
            w_xz[iStep] = JuMP.value(problem.w_xz[iStep])
            w_zy[iStep] = JuMP.value(problem.w_zy[iStep])

            #=
            soc_aux[iStep] = JuMP.value(problem.SOC_aux[iStep])
            p_aux[iStep] = JuMP.value(problem.P_aux[iStep])
            d[iStep] = JuMP.value(problem.d[iStep])       
            deg[iStep] = JuMP.value(problem.deg[iStep])
            d_1[iStep] = JuMP.value(problem.d_1[iStep])
            d_2[iStep] = JuMP.value(problem.d_2[iStep])
            deg_1[iStep] = JuMP.value(problem.deg_1[iStep])
            deg_2[iStep] = JuMP.value(problem.deg_2[iStep])
            u[iStep] = JuMP.value(problem.u[iStep])
            =#

        end

        soc[end] = JuMP.value(problem.soc[end])
        soc_quad[end] = JuMP.value(problem.soc_quad[end])
        x[end] = JuMP.value(problem.x[end])
        y[end] = JuMP.value(problem.y[end])
        z[end] = JuMP.value(problem.z[end])
        #w = JuMP.value(problem.w[iStep])
        w_xx[end] = JuMP.value(problem.w_xx[end])
        w_yy[end] = JuMP.value(problem.w_yy[end])
        w_zz[end] = JuMP.value(problem.w_zz[end])
        w_xy[end] = JuMP.value(problem.w_xy[end])
        w_xz[end] = JuMP.value(problem.w_xz[end])
        w_zy[end] = JuMP.value(problem.w_zy[end])

        for iStage=1:NStages
            soh_final[iStage] = JuMP.value(problem.soh_final[iStage])
            soh_initial[iStage] = JuMP.value(problem.soh_new[iStage])

            deg_stage[iStage] = sum(deg[iStep] for iStep=((iStage-1)*NHoursStage+1):(NHoursStage*iStage))*k
        end

        #
        for iStage=2:(NStages-1)
            revenues_per_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((iStage-1)*NHoursStage+1):(NHoursStage*iStage)) - Battery_price[iStage]*(soh_initial[iStage]-soh_final[iStage-1])
            gain_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((iStage-1)*NHoursStage+1):(NHoursStage*iStage))
            cost_rev[iStage] = Battery_price[iStage]*(soh_initial[iStage]-soh_final[iStage-1])
        end
        #

        revenues_per_stage[1] = sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((1-1)*NHoursStage+1):(NHoursStage*1)) - Battery_price[1]*(soh_initial[1]-min_SOH)
        gain_stage[1]= sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((1-1)*NHoursStage+1):(NHoursStage*1))
        #cost_rev[1] = Battery_price[1]*(soh_initial[1]-min_SOH)-Battery_price[2]*(soh_final[1]-min_SOH)
        cost_rev[1] = Battery_price[1]*(soh_initial[1]-min_SOH)

        
        revenues_per_stage[NStages] = sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((NStages-1)*NHoursStage+1):(NHoursStage*NStages)) + Battery_price[NStages+1]*(soh_final[NStages]-min_SOH)-Battery_price[NStages-1]*(soh_initial[NStages]-soh_final[NStages-1])
        gain_stage[NStages]= sum(Power_prices[iStep]*NHoursStep*max_P*(discharge[iStep]-charge[iStep]) for iStep=((NStages-1)*NHoursStage+1):(NHoursStage*NStages))
        cost_rev[NStages] = -Battery_price[NStages+1]*(soh_final[NStages]-min_SOH) + Battery_price[NStages]*(soh_initial[NStages]-soh_final[NStages-1])
        

    end
    
    println("Collected results")

    return Results(
        objective,
        revenues_per_stage,
        gain_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        deg,
        soc_quad,
        x,
        y,
        z,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
       #= soc_aux,
        p_aux,
        d,
        d_1,
        d_2,
        deg_1,
        deg_2,
        u,
        =#
        soh_final,
        soh_initial,  
    )

end