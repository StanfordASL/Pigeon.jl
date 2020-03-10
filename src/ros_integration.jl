@rosimport std_msgs.msg: ColorRGBA, Header, Float32MultiArray
@rosimport geometry_msgs.msg: Point, TransformStamped, Transform
@rosimport visualization_msgs.msg: Marker
@rosimport safe_traffic_weaving.msg: VehicleTrajectory, path, XYThV
@rosimport auto_messages.msg: from_autobox, to_autobox
rostypegen(@__MODULE__)
import .std_msgs.msg: ColorRGBA, Header, Float32MultiArray
import .geometry_msgs.msg: Point, TransformStamped, Transform
import .visualization_msgs.msg: Marker
import .safe_traffic_weaving.msg: VehicleTrajectory, path, XYThV
import .auto_messages.msg: from_autobox, to_autobox


function fill_2D_transform!(t::TransformStamped, x, y, θ, stamp = RobotOS.now())
    t.header.stamp = stamp
    t.transform.translation.x = x
    t.transform.translation.y = y
    q = Quat(RotZYX(θ, 0, 0))    # tf.transformations[:quaternion_from_euler](0.0, 0.0, θ)
    t.transform.rotation.x = q.x
    t.transform.rotation.y = q.y
    t.transform.rotation.z = q.z
    t.transform.rotation.w = q.w
    t
end

TrajectoryTube(p::path) = TrajectoryTube{Float64}(
    invcumtrapz(p.Ux_des_mps, p.s_m), p.s_m, p.Ux_des_mps, p.Ax_des_mps2,
    p.posE_m, p.posN_m, p.Psi_rad, p.k_1pm, p.grade_rad, 0*p.grade_rad, p.edge_L_m, p.edge_R_m
)
TrajectoryTube(p::VehicleTrajectory) = TrajectoryTube{Float64}(
    p.t, p.s, p.V, p.A, p.E, p.N, p.heading, p.curvature, p.grade, p.bank, p.edge_L, p.edge_R
)

### DEBUG
const show_loginfo = fill(false)

### /des_path and /des_traj
const latest_trajectory = fill(straight_trajectory(30., 5.))
const tracking_mode = fill(:path)
const traj_msg = Array{VehicleTrajectory}(undef)
const from_autobox_msg = Array{from_autobox}(undef)

function nominal_trajectory_callback_coupled(msg, X1CMPC)    # /des_traj
    traj_msg[] = msg
    mpc = X1CMPC
    latest_trajectory[] = TrajectoryTube(msg)
    tracking_mode[] = :traj
    mpc.time_offset = convert(Float64, msg.header.stamp)
    mpc.solved = false
end

function nominal_trajectory_callback_decoupled(msg::path, X1DMPC::TrajectoryTrackingMPC)                 # /des_path
    mpc = X1DMPC
    latest_trajectory[] = TrajectoryTube(msg)
    tracking_mode[] = :path
    mpc.time_offset = NaN
    mpc.solved = false
end


### /to_autobox
const to_autobox_msg = to_autobox()

### /from_autobox
const use_HJI_policy = fill(false)


function from_autobox_callback(msg::from_autobox, to_autobox_pub::Publisher{to_autobox}, X1DMPC::TrajectoryTrackingMPC, X1CMPC::TrajectoryTrackingMPC)
    from_autobox_msg[] = msg
    mpc = (tracking_mode[] == :path ? X1DMPC : X1CMPC)
    mpc.current_state = BicycleState(msg.E_m, msg.N_m, msg.psi_rad, msg.ux_mps, msg.uy_mps, msg.r_radps)
    mpc.current_control = BicycleControl(to_autobox_msg.delta_cmd_rad, to_autobox_msg.fxf_cmd_N, to_autobox_msg.fxr_cmd_N)
    mpc.trajectory = latest_trajectory[]

    show_loginfo[] && RobotOS.logwarn("tracking_mode[] = $tracking_mode")

    if msg.pre_flag == 0
        show_loginfo[] && RobotOS.loginfo("Pigeon MPC: /from_autobox pre_flag == 0, MPC inactive")
        return
    end
    # RobotOS.loginfo("past preflag!")
    if isnan(mpc.time_offset)
        show_loginfo[] && RobotOS.loginfo("Pigeon MPC: time_offset not set, running in path tracking mode")
        _, _, t = path_coordinates(mpc.trajectory, mpc.current_state)
    else
        t = convert(Float64, msg.header.stamp) - mpc.time_offset
        if t < 0 || t > mpc.trajectory.t[end]
            show_loginfo[] && RobotOS.loginfo("Pigeon MPC: current time $t outside trajectory interval [0, $(mpc.trajectory.t[end])]")
            return
        end
    end

    if mpc.current_state.Ux < 1
        show_loginfo[] && RobotOS.loginfo("Pigeon MPC: current speed < 1, pausing MPC while X1 is stopped")
        return
    end

    MPC_steps_missed = msg.header.seq - (mpc.heartbeat + 1)
    if MPC_steps_missed != 0
        RobotOS.logwarn("Pigeon MPC: $(MPC_steps_missed) from_autobox messages lost")
        mpc.heartbeat = msg.header.seq - 1
    end

    t_elapsed = @elapsed begin
        try
            compute_time_steps!(mpc, t)
            compute_linearization_nodes!(mpc)
            update_QP!(mpc)
            solve!(mpc)
        catch err
            RobotOS.logwarn("Pigeon MPC Error: $err\n$(stacktrace(catch_backtrace()))")
        end
    end
    if t_elapsed > 0.01
        RobotOS.logwarn("Pigeon MPC: OSQP took $(1000*t_elapsed) ms at heartbeat $(mpc.heartbeat)")
    else
        show_loginfo[] && RobotOS.loginfo("Pigeon MPC: OSQP took $(1000*t_elapsed) ms at heartbeat $(mpc.heartbeat)")
    end


    mpc.heartbeat += 1

    s, e, _ = path_coordinates(mpc.trajectory, mpc.current_state)

    u_next = get_next_control(mpc)

    to_autobox_msg.header.stamp  = RobotOS.now()
    to_autobox_msg.post_flag     = 1    # TODO: check for OSQP failure
    to_autobox_msg.heartbeat     = mpc.heartbeat
    to_autobox_msg.s_m           = s
    to_autobox_msg.e_m           = e
    to_autobox_msg.delta_cmd_rad = u_next.δ
    to_autobox_msg.fxf_cmd_N     = u_next.Fxf
    to_autobox_msg.fxr_cmd_N     = u_next.Fxr


    
    if isnan(to_autobox_msg.delta_cmd_rad) || isnan(to_autobox_msg.fxf_cmd_N) || isnan(to_autobox_msg.fxr_cmd_N)
        show_loginfo[] && RobotOS.loginfo("Pigeon MPC: OSQP returned NaNs " *
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

### ROS node init
function start_ROS_node(roadway_name="west_paddock", traj_mpc=X1CMPC)
    init_node("pigeon", anonymous=false)

    to_autobox_pub = Publisher{to_autobox}("/to_autobox", queue_size=10)
    RobotOS.loginfo("Publisher to_autobox passed!")
    Subscriber{path}("/des_path", nominal_trajectory_callback_decoupled, (X1DMPC,), queue_size=1)
    Subscriber{VehicleTrajectory}("/des_traj", nominal_trajectory_callback_coupled, (X1CMPC,), queue_size=1)
    Subscriber{from_autobox}("/from_autobox", from_autobox_callback, (to_autobox_pub, X1DMPC, X1CMPC,), queue_size=1)
    RobotOS.loginfo("Subscribers passed!")


    @spawn spin()

end
