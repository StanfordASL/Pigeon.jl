module Pigeon

using LinearAlgebra
using Distributed

using StaticArrays
using DifferentialDynamicsModels
using LinearDynamicsModels
using ForwardDiff
using Interpolations
using OSQP.MathOptInterfaceOSQP
import MathOptInterface
const MOI = MathOptInterface
using Parametron
using RobotOS
using JLD2
# using MAT

import StaticArrays: SUnitRange
import Interpolations: GriddedInterpolation

Parametron.Parameter(A::AbstractArray, model) = Parameter(identity, A, model)

include("math.jl")
include("vehicles.jl")
include("vehicle_dynamics.jl")
include("HJI_computation.jl")
include("trajectories.jl")
include("model_predictive_control.jl")
include("decoupled_lat_long.jl")
include("coupled_lat_long.jl")

const X1DMPC = DecoupledTrajectoryTrackingMPC(X1(), straight_trajectory(30., 5.))
const X1CMPC = CoupledTrajectoryTrackingMPC(X1(), straight_trajectory(30., 5.))
X1DMPC.current_state = BicycleState(0., 0., 0., 5., 0., 0.)
X1DMPC.current_control = BicycleControl(0., 0., 0.)
X1CMPC.current_state = BicycleState(0., 0., 0., 5., 0., 0.)
X1CMPC.current_control = BicycleControl(0., 0., 0.)
X1CMPC.HJI_cache = HJICache(joinpath(@__DIR__, "BicycleCAvoid.jld2"))

include("ros_integration.jl")

function __init__()
    Parametron.initialize!(X1DMPC.model)    # Refresh pointer to OSQP model
    Parametron.initialize!(X1CMPC.model)    # Refresh pointer to OSQP model

    # dry-runs to force compilation
    compute_time_steps!(X1DMPC, 0.)
    compute_linearization_nodes!(X1DMPC)
    update_QP!(X1DMPC)
    solve!(X1DMPC)

    compute_time_steps!(X1CMPC, 0.)
    compute_linearization_nodes!(X1CMPC)
    update_QP!(X1CMPC)
    solve!(X1CMPC)
end

function load_test_paths()
    include(joinpath(@__DIR__, "..", "test", "load_test_paths.jl"))
end

end # module
