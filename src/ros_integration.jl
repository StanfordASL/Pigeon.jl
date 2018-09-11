@rosimport osprey.msg: path
@rosimport asl_prototyping.msg: VehicleTrajectory, XYThV
@rosimport auto_messages.msg: from_autobox, to_autobox
rostypegen(@__MODULE__)
import .osprey.msg: path
import .asl_prototyping.msg: VehicleTrajectory, XYThV
import .auto_messages.msg: from_autobox, to_autobox

TrajectoryTube(p::path) = TrajectoryTube{Float64}(
    invcumtrapz(p.Ux_des_mps, p.s_m), p.s_m, p.Ux_des_mps, p.Ax_des_mps2,
    p.posE_m, p.posN_m, p.Psi_rad, p.k_1pm, p.grade_rad, 0*p.grade_rad, p.edge_L_m, p.edge_R_m
)
TrajectoryTube(p::VehicleTrajectory) = TrajectoryTube{Float64}(
    p.t, p.s, p.V, p.A, p.E, p.N, p.heading, p.curvature, p.grade, p.bank, p.edge_L, p.edge_R
)

### /des_path and /des_traj
const latest_trajectory = fill(straight_trajectory(30., 5.))
const tracking_mode = fill(:path)
function nominal_trajectory_callback(msg::path, mpc=X1DMPC)                 # /des_traj
    latest_trajectory[] = TrajectoryTube(msg)
    tracking_mode[] = :path
    mpc.time_offset = NaN
    mpc.solved = false
end
function nominal_trajectory_callback(msg::VehicleTrajectory, mpc=X1CMPC)    # /des_path
    latest_trajectory[] = TrajectoryTube(msg)
    tracking_mode[] = :traj
    mpc.time_offset = convert(Float64, msg.header.stamp)
    mpc.solved = false
end

### /to_autobox
const to_autobox_msg = to_autobox()

### /from_autobox
const use_HJI_policy = fill(false)
function from_autobox_callback(msg::from_autobox, to_autobox_pub, path_mpc=X1DMPC, traj_mpc=X1CMPC)
    mpc = (tracking_mode[] == :path ? path_mpc : traj_mpc)
    mpc.current_state = BicycleState(msg.E_m, msg.N_m, msg.psi_rad, msg.ux_mps, msg.uy_mps, msg.r_radps)
    # mpc.current_control = BicycleControl(msg.delta_rad, msg.fxf_N, msg.fxr_N)
    mpc.current_control = BicycleControl(to_autobox_msg.delta_cmd_rad, to_autobox_msg.fxf_cmd_N, to_autobox_msg.fxr_cmd_N)
    mpc.trajectory = latest_trajectory[]
    if msg.pre_flag == 0
        RobotOS.loginfo("Pigeon MPC: /from_autobox pre_flag == 0, MPC inactive")
        return
    end
    if isnan(mpc.time_offset)
        RobotOS.loginfo("Pigeon MPC: time_offset not set, running in path tracking mode")
        _, _, t = path_coordinates(mpc.trajectory, mpc.current_state)
    else
        t = convert(Float64, msg.header.stamp) - mpc.time_offset
        if t < 0 || t > mpc.trajectory.t[end]
            RobotOS.loginfo("Pigeon MPC: current time $t outside trajectory interval [0, $(mpc.trajectory.t[end])]")
            return
        end
    end
    if mpc.current_state.Ux < 1
        RobotOS.loginfo("Pigeon MPC: current speed < 1, pausing MPC while X1 is stopped")
        return
    end
    MPC_steps_missed = msg.header.seq - (mpc.heartbeat + 1)
    if MPC_steps_missed != 0
        RobotOS.logwarn("Pigeon MPC: $(MPC_steps_missed) from_autobox messages lost")
        mpc.heartbeat = msg.header.seq - 1
    end

    local relative_state, V, ∇V
    t_elapsed = @elapsed begin
        try
            if tracking_mode[] == :traj
                relative_state = HJIRelativeState(mpc.current_state, mpc.other_car_state)
                V, ∇V = mpc.HJI_cache[relative_state]
                RobotOS.loginfo("Pigeon MPC: HJI value function = $V")
                # RobotOS.loginfo("Pigeon MPC: HJI relative state = $(relative_state)")
            end
            compute_time_steps!(mpc, t)
            compute_linearization_nodes!(mpc)
            update_QP!(mpc)
            solve!(mpc)
        catch err
            RobotOS.logwarn("Pigeon MPC Error: $err\n$(stacktrace(catch_backtrace()))")
        end
    end

    RobotOS.loginfo("Pigeon MPC: OSQP took $(1000*t_elapsed) ms at heartbeat $(mpc.heartbeat)")
    # RobotOS.loginfo("Pigeon MPC: $(mpc.model.optimizer.results.info)")
    # RobotOS.loginfo("deltas: $(value.(mpc.model, mpc.variables.δ))")
    mpc.heartbeat += 1

    s, e, _ = path_coordinates(mpc.trajectory, mpc.current_state)
    if tracking_mode[] == :traj && use_HJI_policy[] && V <= mpc.HJI_ϵ
        u_next = BicycleControl(mpc.dynamics.longitudinal_params, optimal_control(mpc.dynamics, relative_state, ∇V))
        RobotOS.loginfo("Pigeon MPC: HJI stepping in to save the day (with a hammer)")
    else
        if tracking_mode[] == :traj && V <= mpc.HJI_ϵ
            RobotOS.loginfo("Pigeon MPC: HJI stepping in to save the day (with a feather)")
        end
        u_next = get_next_control(mpc)
    end
    to_autobox_msg.header.stamp  = RobotOS.now()
    to_autobox_msg.post_flag     = 1    # TODO: check for OSQP failure
    to_autobox_msg.heartbeat     = mpc.heartbeat
    to_autobox_msg.s_m           = s
    to_autobox_msg.e_m           = e
    to_autobox_msg.delta_cmd_rad = u_next.δ
    to_autobox_msg.fxf_cmd_N     = u_next.Fxf
    to_autobox_msg.fxr_cmd_N     = u_next.Fxr
    if isnan(to_autobox_msg.delta_cmd_rad) || isnan(to_autobox_msg.fxf_cmd_N) || isnan(to_autobox_msg.fxr_cmd_N)
        RobotOS.loginfo("Pigeon MPC: OSQP returned NaNs " *
                        "($(to_autobox_msg.delta_cmd_rad), $(to_autobox_msg.fxf_cmd_N), $(to_autobox_msg.fxr_cmd_N))" *
                        "; using previous control "*
                        "($(mpc.current_control.δ), $(mpc.current_control.Fxf), $(mpc.current_control.Fxr))")
        to_autobox_msg.delta_cmd_rad = mpc.current_control.δ
        to_autobox_msg.fxf_cmd_N     = mpc.current_control.Fxf
        to_autobox_msg.fxr_cmd_N     = mpc.current_control.Fxr
        publish(to_autobox_pub, to_autobox_msg)
        to_autobox_msg.delta_cmd_rad = 0    # ensures that if we get 2 NaNs in a row, we return to 0 control
        to_autobox_msg.fxf_cmd_N     = 0
        to_autobox_msg.fxr_cmd_N     = 0
        Parametron.initialize!(mpc.model)   # reinitialize model (otherwise NaNs might stick)
        mpc.solved = false    # any warm starting should ignore the previous solution
    else
        publish(to_autobox_pub, to_autobox_msg)
    end
end

function other_car_callback(msg::XYThV, mpc=X1CMPC)
    mpc.other_car_state = SimpleCarState(msg.x, msg.y, msg.th - pi/2, msg.v)
end

### ROS node init
function start_ROS_node()
    init_node("pigeon", anonymous=false)
    to_autobox_pub = Publisher{to_autobox}("/to_autobox", queue_size=10)
    Subscriber{path}("/des_path", nominal_trajectory_callback, queue_size=1)
    Subscriber{VehicleTrajectory}("/des_traj", nominal_trajectory_callback, queue_size=1)
    Subscriber{from_autobox}("/from_autobox", from_autobox_callback, (to_autobox_pub,), queue_size=1)
    Subscriber{XYThV}("/xbox_car/xythv", other_car_callback, queue_size=1)    # TODO: abstract with a republisher
    @spawn spin()
end
