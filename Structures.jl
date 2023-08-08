# STRUCTURES USED IN THE PROBLEM

# Input data
#-----------------------------------------------

# Input parameters 
@with_kw struct InputParam{F<:Float64,I<:Int}
    NYears::F                                     # Number of years
    NMonths::I
    NStages::I                                    # Number of stages of N months in the problem FORMULATION-- calcolato come NYears/NMonths*12
    NHoursStep::F                                 # Number of hours in each time step 
    NHoursStage::I                                # Number of hours in each Stage (3-4-6 months)
    NSteps::I                                     # Number of steps in the NYeras --> NYears*8760/NHoursStep
    Big::F                                        # A big number
    conv::F                                       # A small number for degradation convergence
    disc::I                                       # Discretization points
end

# Battery's characteristics
@with_kw struct BatteryParam{F<:Float64,I<:Int}
    min_SOC::F                                     # Battery's minimum energy storage capacity
    max_SOC::F                                     # Batter's maximum energy storage capacity
    Eff_charge::F                                  # Battery's efficiency for charging operation
    Eff_discharge::F                               # Battery's efficiency for discharging operation
    min_P::F
    max_P::F
    max_SOH::F                                     # Maximum SOH that can be achieved because of volume issues
    min_SOH::F                                     # Minimum SOH to be respected by contract
    Nfull::I                                       # Maximum number of full cycles for DoD=100%                                        
end
  
# solver parameters
@with_kw struct SolverParam{F<:Float64,I<:Int}
    MIPGap::F 
    MIPFocus::I
    Method::F
    Cuts::F
    Heuristics::F
end

#= envelopes
struct envParam{F<:Float64,I<:Int, S<:String}
    n_envelopes::I
    low_x::Any
    up_x::Any
end
=#
  
# Indirizzi cartelle
@with_kw struct caseData{S<:String}
    DataPath::S
    InputPath::S
    ResultPath::S
    CaseName::S
end

# runMode Parameters
@with_kw mutable struct runModeParam{B<:Bool}

    # Solver settings
    solveMIP::B     #If using SOS2

    batterySystemFromFile::B 

    #runMode self defined reading of input 
    setInputParameters::B             #from .in file
 
    excel_savings::B 

end

# Optimization problem
struct BuildStageProblem
    M::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    w_uu::Any
    w_xu::Any
    w_yu::Any
    w_zu::Any
    soh_final::Any
    soh_new::Any
end

struct Results
    objective::Any
    revenues_per_stage::Any
    gain_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    charge::Any
    discharge::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_uu::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    w_xu::Any
    w_yu::Any
    w_zu::Any
    soh_final::Any
    soh_initial::Any
end
