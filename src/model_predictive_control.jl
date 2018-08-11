struct MPCTimeSteps{T}
    ts::Vector{T}
    dt::Vector{T}
    N_short::Int
    N_long::Int
    dt_short::T
    dt_long::T
    use_correction_step::Bool
    prev_ts::Vector{T}
end
function MPCTimeSteps(N_short::Int, N_long::Int, dt_short::T, dt_long::T, use_correction_step::Bool) where {T}
    N = 1 + N_short + N_long
    ts = T.(1:N)    # initialized so that dt's are nonzero
    dt = diff(ts)
    MPCTimeSteps(ts, dt, N_short, N_long, dt_short, dt_long, use_correction_step, ts)
end
function compute_time_steps!(TS::MPCTimeSteps, t0)
    N_short, N_long, dt_short, dt_long, use_correction_step = TS.N_short, TS.N_long, TS.dt_short, TS.dt_long, TS.use_correction_step

    TS.prev_ts .= TS.ts
    t0_long = t0 + N_short*dt_short
    if use_correction_step
        t0_long = dt_long*ceil((t0_long+dt_short)/dt_long - 1)    # dt_correction ∈ [dt_short, dt_long+dt_short)
    end
    TS.ts[1:N_short+1] .= t0 .+ dt_short*(0:N_short)
    TS.ts[N_short+2:end] .= t0_long .+ dt_long*(1:N_long)
    for i in 1:N_short+N_long    # long-winded diff!
        TS.dt[i] = TS.ts[i+1] - TS.ts[i]
    end
end

mutable struct TrajectoryTrackingMPC{T,C,Q,U,P,QPP,QPV}
    vehicle::Dict{Symbol,T}
    trajectory::TrajectoryTube{T}
    dynamics::VehicleModel{T}
    control_params::C

    current_state::BicycleState{T}
    current_control::BicycleControl{T}
    heartbeat::Int
    time_offset::Float64

    time_steps::MPCTimeSteps{T}

    qs::Vector{Q}
    us::Vector{U}
    ps::Vector{P}

    tracking_dynamics::VehicleModel{T}
    model::Model{T,OSQPOptimizer}
    variables::QPV
    parameters::QPP
    solved::Bool
end

compute_time_steps!(mpc::TrajectoryTrackingMPC, t0) = compute_time_steps!(mpc.time_steps, t0)

compute_linearization_nodes!(mpc::TrajectoryTrackingMPC) = compute_linearization_nodes!(mpc, mpc.qs, mpc.us, mpc.ps)

update_QP!(mpc::TrajectoryTrackingMPC) = update_QP!(mpc, mpc.parameters)

Parametron.solve!(mpc::TrajectoryTrackingMPC) = (solve!(mpc.model); mpc.solved = true)

get_next_control(mpc::TrajectoryTrackingMPC) = get_next_control(mpc, mpc.variables)

function simulate(mpc, q0, u0, dt=0.01)
    mpc.current_state = q0
    mpc.current_control = u0
    qs = []
    xs = []
    us = []
    ps = []
    for t in 0:dt:mpc.trajectory.t[end]
        push!(qs, mpc.current_state)
        push!(us, mpc.current_control)
        compute_time_steps!(mpc, t)
        compute_linearization_nodes!(mpc)
        update_QP!(mpc)
        solve!(mpc)
        mpc.current_state = propagate(mpc.dynamics, mpc.current_state, StepControl(dt, BicycleControl2(mpc.current_control)))
        mpc.current_control = get_next_control(mpc)
        push!(xs, mpc.qs[1])
        push!(ps, mpc.ps[1])
    end
    qs, xs, us, ps
end
