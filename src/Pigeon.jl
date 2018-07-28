__precompile__()

module Pigeon

using LinearAlgebra
using Distributed
using StaticArrays
using DifferentialDynamicsModels
using Interpolations
using ForwardDiff
using OSQP.MathOptInterfaceOSQP
import MathOptInterface
const MOI = MathOptInterface
using SimpleQP
using RobotOS

macro ST(expr)    # ensures that subtypes of StaticArrays.FieldVector maintain their type upon manipulation
    full_name = expr.args[2].args[1]
    short_name = full_name isa Symbol ? full_name : full_name.args[1]
    N = expr.args[2].args[2].args[2]
    T = expr.args[2].args[2].args[3]
    quote
        $(esc(expr))
        $(esc(:(StaticArrays.similar_type(::Type{<:$short_name}, ::Type{$T}, s::Size{($N,)}) where {$T} = $full_name)))
    end
end

include("math.jl")
include("trajectories.jl")
include("vehicles.jl")
include("vehicle_dynamics.jl")
include("qp.jl")
include("model_predictive_control.jl")

const X1MPC = MPC(X1(), straight_trajectory(30., 5.))
### dry-run to force compilation
X1MPC.current_state = BicycleState(0., 0., 0., 5., 0., 0.)
X1MPC.current_control = BicycleControl(0., 0., 0.)
MPC_time_steps!(X1MPC, 0.)
compute_linearization_nodes!(X1MPC)
update_QP!(X1MPC)
solve!(X1MPC.model)

include("ros_integration.jl")

function __init__()
    @rosimport osprey.msg: path                             # PyObject initialization
    @rosimport asl_prototyping.msg: VehicleTrajectory
    @rosimport auto_bridge.msg: from_autobox, to_autobox
    SimpleQP.initialize!(X1MPC.model)    # Refresh pointer to OSQP model

    # more compilation
    MPC_time_steps!(X1MPC, 0.)
    compute_linearization_nodes!(X1MPC)
    update_QP!(X1MPC)
    solve!(X1MPC.model)
end

function load_test_paths()
    include(joinpath(@__DIR__, "..", "test", "load_test_paths.jl"))
end

end # module
