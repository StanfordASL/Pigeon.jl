__precompile__()

module Pigeon

using StaticArrays
using DifferentialDynamicsModels
using Interpolations
using ForwardDiff
using OSQP
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

function __init__()

@eval Pigeon begin

SimpleQP.initialize!(X1MPC.model)    # Refresh pointer to OSQP model

@rosimport osprey.msg: path
@rosimport asl_prototyping.msg: VehicleTrajectory
@rosimport auto_bridge.msg: from_autobox, to_autobox
rostypegen()
import osprey.msg: path
import asl_prototyping.msg: VehicleTrajectory
import auto_bridge.msg: from_autobox, to_autobox

TrajectoryTube(p::path) = TrajectoryTube{Float64}(
    invcumtrapz(p.Ux_des_mps, p.s_m), p.s_m, p.Ux_des_mps, p.Ax_des_mps2,
    p.posE_m, p.posN_m, p.Psi_rad, p.k_1pm, p.grade_rad, 0*p.grade_rad, p.edge_L_m, p.edge_R_m
)
TrajectoryTube(p::VehicleTrajectory) = TrajectoryTube{Float64}(
    p.t, p.s, p.V, p.A, p.E, p.N, p.heading, p.curvature, p.grade, p.bank, p.edge_L, p.edge_R
)

### /des_path
const latest_trajectory = fill(straight_trajectory(30., 5.))
function nominal_trajectory_callback(msg::path, mpc=X1MPC)
    latest_trajectory[] = TrajectoryTube(msg)
    if isnan(mpc.time_offset)
        mpc.time_offset = Float64(msg.header.stamp)
    end
end
function nominal_trajectory_callback(msg::VehicleTrajectory, mpc=X1MPC)
    latest_trajectory[] = TrajectoryTube(msg)
    if isnan(mpc.time_offset)
        mpc.time_offset = Float64(msg.header.stamp)
    end
end

### /to_autobox
const to_autobox_msg = to_autobox()
const to_autobox_pub = Publisher{to_autobox}("/to_autobox", queue_size=10)

### /from_autobox
const latest_from_autobox = fill(from_autobox())
function from_autobox_callback(msg::from_autobox, mpc=X1MPC)
    mpc.current_state = BicycleState(msg.E_m, msg.N_m, msg.psi_rad, msg.ux_mps, msg.uy_mps, msg.r_radps)
    # mpc.current_control = BicycleControl(msg.delta_rad, msg.fxf_N, msg.fxr_N)
    mpc.current_control = BicycleControl(to_autobox_msg.delta_cmd_rad, to_autobox_msg.fxf_cmd_N, to_autobox_msg.fxr_cmd_N)
    mpc.trajectory = latest_trajectory[]
    if isnan(mpc.time_offset)
        RobotOS.loginfo("Pigeon MPC: time_offset not set")
        return
    end
    t = Float64(msg.header.stamp) - mpc.time_offset
    if t < 0 || t > mpc.trajectory.t[end]
        RobotOS.loginfo("Pigeon MPC: current time $t outside range [0, $(mpc.trajectory.t[end])]")
        return
    end
    MPC_steps_missed = msg.header.seq - (mpc.heartbeat + 1)
    if MPC_steps_missed != 0
        RobotOS.logwarn("Pigeon MPC: $(MPC_steps_missed) from_autobox messages lost")
        mpc.heartbeat = msg.header.seq - 1
    end

    tic()
    try
        MPC_time_steps!(mpc, t)
        compute_linearization_nodes!(mpc)
        update_QP!(mpc)
        solve!(mpc.model)
    catch err
        RobotOS.logwarn("OSQP Error: $err")
    end
    t_elapsed = toq()

    RobotOS.loginfo("Pigeon MPC: OSQP took $(1000*t_elapsed) ms at heartbeat $(mpc.heartbeat)")
    RobotOS.loginfo("Pigeon MPC: $(mpc.model.optimizer.results.info)")
    RobotOS.loginfo("deltas: $(value.(mpc.model, mpc.variables.δ))")
    mpc.heartbeat += 1

    s, e, _ = path_coordinates(mpc.trajectory, SVector(mpc.current_state.E, mpc.current_state.N))
    to_autobox_msg.header.stamp  = RobotOS.now()
    to_autobox_msg.post_flag     = 1    # TODO: check for OSQP failure
    to_autobox_msg.heartbeat     = mpc.heartbeat
    to_autobox_msg.s_m           = s
    to_autobox_msg.e_m           = e
    to_autobox_msg.delta_cmd_rad = value(mpc.model, mpc.variables.δ[2])
    to_autobox_msg.fxf_cmd_N     = mpc.us[2].Fxf
    to_autobox_msg.fxr_cmd_N     = mpc.us[2].Fxr
    publish(to_autobox_pub, to_autobox_msg)
end

### ROS node init
function start_ROS_node()
    init_node("pigeon", anonymous=false)
    Subscriber{path}("/des_path", nominal_trajectory_callback, queue_size=1)
    Subscriber{VehicleTrajectory}("/des_traj", nominal_trajectory_callback, queue_size=1)
    Subscriber{from_autobox}("/from_autobox", from_autobox_callback, queue_size=1)
    @spawn spin()
end

end # eval block

end # __init__

end # module
