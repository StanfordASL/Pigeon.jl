mutable struct MPC{T}
    vehicle::Dict{Symbol,T}
    trajectory::TrajectoryTube{T}
    bicycle_model::BicycleModel{T}
    control_params::ControlParams{T}

    current_state::BicycleState{T}
    current_control::BicycleControl{T}
    heartbeat::Int
    time_offset::Float64

    N_short::Int
    N_long::Int
    dt_short::T
    dt_long::T
    use_correction_step::Bool

    ts::Vector{T}
    dt::Vector{T}
    qs::Vector{TrackingBicycleState{T}}
    us::Vector{BicycleControl{T}}
    ps::Vector{TrackingBicycleParams{T}}

    model::Model{T,OSQPOptimizer}
    variables::QP_Variables
    parameters::QP_Parameters{T}

    function MPC{T}(vehicle::Dict{Symbol,T}, trajectory::TrajectoryTube{T}, control_params::ControlParams{T},
                    N_short::Int, N_long::Int, dt_short::T, dt_long::T, use_correction_step::Bool) where {T}
        bicycle_model = BicycleModel(vehicle)
        current_state = zeros(BicycleState{T})
        current_control = zeros(BicycleControl{T})
        N = 1 + N_short + N_long
        ts = T.(1:N) # zeros(T, N) (not zeros so that construct_QP below works)
        dt = diff(ts)
        qs = rand(TrackingBicycleState{T}, N) # not zeros so that construct_QP below works
        us = rand(BicycleControl{T}, N)
        ps = rand(TrackingBicycleParams{T}, N)
        mpc = new(vehicle, trajectory, bicycle_model, control_params,
                  current_state, current_control, 0, NaN,
                  N_short, N_long, dt_short, dt_long, use_correction_step,
                  ts, dt, qs, us, ps)
        mpc.model, mpc.variables, mpc.parameters = construct_QP(mpc)
        mpc
    end
end
function MPC(vehicle::Dict{Symbol,T}, trajectory::TrajectoryTube{T}; control_params=ControlParams(vehicle),
             N_short=10, N_long=20, dt_short=T(0.01), dt_long=T(0.2), use_correction_step=true) where {T}
    MPC{T}(vehicle, trajectory, control_params, N_short, N_long, dt_short, dt_long, use_correction_step)
end

### Time Steps
function MPC_time_steps!(ts::Vector, dt::Vector, t0, N_short=10, N_long=20, dt_short=0.01, dt_long=0.2, use_correction_step=true)
    t0_long = t0 + N_short*dt_short
    if use_correction_step
        t0_long = dt_long*ceil((t0_long+dt_short)/dt_long - 1)    # dt_correction ∈ [dt_short, dt_long+dt_short)
    end
    ts[1:N_short+1] .= t0 + dt_short*(0:N_short)
    ts[N_short+2:end] .= t0_long + dt_long*(1:N_long)
    for i in 1:N_short+N_long    # long-winded diff!
        dt[i] = ts[i+1] - ts[i]
    end
    ts
end
MPC_time_steps!(mpc::MPC, t0) = MPC_time_steps!(mpc.ts, mpc.dt, t0, mpc.N_short, mpc.N_long, mpc.dt_short, mpc.dt_long, mpc.use_correction_step)
MPC_time_steps(t0, N_short=10, N_long=20, dt_short=0.01, dt_long=0.2, use_correction_step=true) =
    MPC_time_steps!(zeros(typeof(t0), 1+N_short+N_long), zeros(typeof(t0), N_short+N_long), t0, N_short, N_long, dt_short, dt_long, use_correction_step)

function compute_linearization_nodes!(mpc::MPC{T}) where {T}
    traj = mpc.trajectory
    B = mpc.bicycle_model
    U = mpc.control_params
    q0 = mpc.current_state
    u0 = mpc.current_control
    N_short, N_long = mpc.N_short, mpc.N_long
    ts, dt = mpc.ts, mpc.dt

    s, e0, t0 = path_coordinates(traj, SVector(q0.E, q0.N))

    V = hypot(q0.Ux, q0.Uy)
    β0 = atan2(q0.Uy, q0.Ux)
    r0 = q0.r
    δ0 = u0.δ
    Fyf0, Fyr0 = lateral_tire_forces(B, q0, u0)
    for i in 1:N_short+N_long+1
        τ = (i == N_short+N_long+1 ? dt[i-1] : dt[i])
        tj = traj[s]
        κ = tj.κ
        A_des = tj.A + U.k_V*(tj.V - V)/τ + U.k_s*(traj(ts[i]).s - s)/τ/τ   # alternatively, - U.k_t*(tj.t - ts[i])/τ/τ
        if i <= N_short+1
            q = TrackingBicycleState(q0.Uy, q0.r, adiff(q0.ψ, tj.ψ), e0)    # to match paper should be (..., 0.0, 0.0)
            _, u, p, A = steady_state_estimates(B, U, κ, A_des, V, 1, r0, β0, δ0, Fyf0)
        else
            q, u, p, A = steady_state_estimates(B, U, κ, A_des, V)
        end
        mpc.qs[i] = q
        mpc.us[i] = u
        mpc.ps[i] = p
        i == N_short+N_long+1 && break
        V = V + A*τ
        s = s + V*τ + A*τ*τ
    end
end

function update_QP!(mpc)
    BM = mpc.bicycle_model
    U = mpc.control_params
    N_short, N_long = mpc.N_short, mpc.N_long
    dt = mpc.dt
    qs, us, ps = mpc.qs, mpc.us, mpc.ps
    curr_δ = mpc.current_control.δ
    QPP = mpc.parameters

    QPP.Q_Δψ().diag .= U.Q_Δψ .* dt
    QPP.Q_e().diag  .= U.Q_e .* dt
    QPP.R_δ().diag  .= U.R_δ .* dt
    QPP.R_Δδ().diag .= U.R_Δδ ./ dt
    QPP.W_β()    .= U.W_β .* dt
    QPP.W_r()    .= U.W_r .* dt
    QPP.curr_q() .= qs[1]
    QPP.curr_δ() .= curr_δ

    for t in 1:N_short
        At, Bt, ct = ZOH(BM, qs[t], us[t], ps[t], dt[t])
        QPP.A[t]() .= At
        QPP.B[t]() .= Bt
        QPP.c[t]() .= ct
    end
    for t in N_short+1:N_short+N_long
        At, B0t, B1t, ct = FOH(BM, qs[t], us[t], ps[t], us[t+1], ps[t+1], dt[t])
        QPP.A[t]() .= At
        QPP.B0[t-N_short]() .= B0t
        QPP.B1[t-N_short]() .= B1t
        QPP.c[t]() .= ct
    end
    δ_hwmax, δ̇_max = U.δ_hwmax, U.δ̇_max
    for t in 1:N_short+N_long
        δ_mint, δ_maxt, Ht, Gt = stable_limits(BM, qs[t+1], us[t+1], ps[t+1])
        QPP.H[t]() .= Ht
        QPP.G[t]() .= Gt
        QPP.δ_min[t]() .= max(δ_mint, -δ_hwmax)
        QPP.δ_max[t]() .= min(δ_maxt, δ_hwmax)
        QPP.Δδ_min[t]() .= -δ̇_max*dt[t]
        QPP.Δδ_max[t]() .= δ̇_max*dt[t]
    end
end

function simulate(mpc, q0, u0, dt=0.01)
    mpc.current_state = q0
    mpc.current_control = u0
    qs = []
    us = []
    for t in 0:dt:mpc.trajectory.t[end]
        push!(qs, mpc.current_state)
        push!(us, mpc.current_control)
        MPC_time_steps!(mpc, t)
        compute_linearization_nodes!(mpc)
        update_QP!(mpc)
        solve!(mpc.model)
        mpc.current_state = propagate(mpc.bicycle_model, mpc.current_state, StepControl(dt, mpc.current_control))
        mpc.current_control = BicycleControl(value(mpc.model, mpc.variables.δ[2]), mpc.us[2][2], mpc.us[2][3])
    end
    qs, us
end
