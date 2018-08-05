module Pigeon

using LinearAlgebra
using Distributed

using StaticArrays
using DifferentialDynamicsModels
using LinearDynamicsModels
using Interpolations
using OSQP.MathOptInterfaceOSQP
import MathOptInterface
const MOI = MathOptInterface
using SimpleQP
using RobotOS

SimpleQP.Parameter(A::AbstractArray, model) = Parameter(identity, A, model)

include("math.jl")
include("vehicles.jl")
include("vehicle_dynamics.jl")
include("trajectories.jl")
include("model_predictive_control.jl")
include("decoupled_lat_long.jl")
# include("coupled_lat_long.jl")

const X1MPC = DecoupledTrajectoryTrackingMPC(X1(), straight_trajectory(30., 5.))
### dry-run to force compilation
X1MPC.current_state = BicycleState(0., 0., 0., 5., 0., 0.)
X1MPC.current_control = BicycleControl(0., 0., 0.)
compute_time_steps!(X1MPC, 0.)
compute_linearization_nodes!(X1MPC)
update_QP!(X1MPC)
solve!(X1MPC)

include("ros_integration.jl")

function __init__()
    SimpleQP.initialize!(X1MPC.model)    # Refresh pointer to OSQP model

    # more compilation
    compute_time_steps!(X1MPC, 0.)
    compute_linearization_nodes!(X1MPC)
    update_QP!(X1MPC)
    solve!(X1MPC)
end

function load_test_paths()
    include(joinpath(@__DIR__, "..", "test", "load_test_paths.jl"))
end

end # module
