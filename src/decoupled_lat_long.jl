struct DecoupledControlParams{T}
    V_min::T
    V_max::T

    k_V::T
    k_s::T

    δ̇_max::T

    Q_Δψ::T
    Q_e::T
    W_β::T
    W_r::T

    R_δ::T
    R_Δδ::T
end
function DecoupledControlParams(;V_min=2.0,
                                 V_max=15.0,
                                 k_V=10/4/100,
                                 k_s=10/4/10000,
                                 δ̇_max=0.344,
                                 Q_Δψ=1/(10*π/180)^2,
                                 Q_e=1.0,
                                 W_β=50/(10*π/180),
                                 W_r=50.0,
                                 R_δ=0.0,
                                 R_Δδ=0.01/(10*π/180)^2)
    DecoupledControlParams(V_min, V_max, k_V, k_s, δ̇_max, Q_Δψ, Q_e, W_β, W_r, R_δ, R_Δδ)
end

function DecoupledTrajectoryTrackingMPC(vehicle::Dict{Symbol,T}, trajectory::TrajectoryTube{T}; control_params=DecoupledControlParams(),
                                        N_short=10, N_long=20, dt_short=T(0.01), dt_long=T(0.2), use_correction_step=true) where {T}
    dynamics = VehicleModel(vehicle)
    current_state = zeros(BicycleState{T})
    current_control = zeros(BicycleControl{T})
    time_steps = MPCTimeSteps(N_short, N_long, dt_short, dt_long, use_correction_step)

    N = 1 + N_short + N_long
    qs = rand(LateralTrackingBicycleState{T}, N)    # not zeros so that construct_lateral_tracking_QP below works
    us = rand(BicycleControl2{T}, N)
    ps = rand(LateralTrackingBicycleParams{T}, N)
    tracking_dynamics = VehicleModel(vehicle, LateralTrackingBicycleModel(vehicle))
    model, variables, parameters = construct_lateral_tracking_QP(tracking_dynamics, control_params, time_steps, qs, us, ps)
    mpc = TrajectoryTrackingMPC(vehicle, trajectory, dynamics, control_params,
                                current_state, current_control, 0, NaN,
                                time_steps,
                                qs, us, ps,
                                tracking_dynamics, model, variables, parameters, false)
    mpc
end

function compute_linearization_nodes!(mpc::TrajectoryTrackingMPC{T},
                                      qs::Vector{LateralTrackingBicycleState{T}},
                                      us::Vector{BicycleControl2{T}},
                                      ps::Vector{LateralTrackingBicycleParams{T}}) where {T}
    traj = mpc.trajectory
    X  = mpc.tracking_dynamics
    U  = mpc.control_params
    q0 = mpc.current_state
    u0 = mpc.current_control
    TS = mpc.time_steps
    N_short, N_long = TS.N_short, TS.N_long
    ts, dt = TS.ts, TS.dt

    s, e0, t0 = path_coordinates(traj, q0)

    V = hypot(q0.Ux, q0.Uy)
    β0 = atan(q0.Uy, q0.Ux)
    r0 = q0.r
    δ0 = u0.δ
    Fyf0, Fyr0 = lateral_tire_forces(X.bicycle_model, q0, u0)
    for i in 1:N_short+N_long+1
        τ = (i == N_short+N_long+1 ? dt[i-1] : dt[i])
        tj = traj[s]
        κ = tj.κ
        A_des = tj.A + U.k_V*(tj.V - V)/τ + (isnan(mpc.time_offset) ? 0 : U.k_s*(traj(ts[i]).s - s)/τ/τ)   # alternatively, - U.k_t*(tj.t - ts[i])/τ/τ
        A_des = min(max(A_des, (U.V_min - V)/τ), (U.V_max - V)/τ)
        if i == 1
            q = LateralTrackingBicycleState(q0.Uy, q0.r, adiff(q0.ψ, tj.ψ), e0)    # to match paper should be (..., 0.0, 0.0)
            u = BicycleControl2(u0)
            p = LateralTrackingBicycleParams(q0.Ux, κ, T(0), T(0))
            q̇ = mpc.dynamics(q0, BicycleControl2(u0), LocalRoadGeometry(tj))
            A = (q̇[4] - q0.r*q0.Uy)*cos(β0) + (q̇[5] + q0.r*q0.Ux)*sin(β0)
        elseif i <= N_short+1
            q = LateralTrackingBicycleState(q0.Uy, q0.r, adiff(q0.ψ, tj.ψ), e0)    # to match paper should be (..., 0.0, 0.0)
            est = steady_state_estimates(X, V, A_des, κ, num_iters=1, r=r0, β0=β0, δ0=δ0, Fyf0=Fyf0)
            u = BicycleControl2(est.δ, est.Fxf + est.Fxr)
            p = LateralTrackingBicycleParams(est.Ux, κ, T(0), T(0))
            A = est.A
        else
            est = steady_state_estimates(X, V, A_des, κ)
            q = LateralTrackingBicycleState(est.Uy, est.r, -est.β, T(0))
            u = BicycleControl2(est.δ, est.Fxf + est.Fxr)
            p = LateralTrackingBicycleParams(est.Ux, κ, T(0), T(0))
            A = est.A
        end
        qs[i] = q
        us[i] = u
        ps[i] = p
        i == N_short+N_long+1 && break
        V = V + A*τ
        s = s + V*τ + A*τ*τ
    end
end

struct LateralTrackingQPParams{T}
    Q_Δψ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    Q_e   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_δ   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_Δδ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    W_β   ::Parameter{Vector{T},typeof(identity),true}
    W_r   ::Parameter{Vector{T},typeof(identity),true}
    q_curr::Parameter{Vector{T},typeof(identity),true}
    δ_curr::Parameter{Vector{T},typeof(identity),true}
    A     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B0    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    Bf    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    c     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    H     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    G     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_min ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_max ::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_min::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_max::Vector{Parameter{Vector{T},typeof(identity),true}}
end

struct LateralTrackingQPVariables
    q::Matrix{Variable}
    δ::Matrix{Variable}
    σ::Matrix{Variable}
end

function construct_lateral_tracking_QP(dynamics::VehicleModel{T}, control_params, time_steps, qs, us, ps) where {T}
    N_short, N_long, dt = time_steps.N_short, time_steps.N_long, time_steps.dt

    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.WarmStart(), true)
    m = Model{T}(optimizer)

    Q_Δψ   = Parameter(Diagonal(control_params.Q_Δψ .* dt), m)
    Q_e    = Parameter(Diagonal(control_params.Q_e  .* dt), m)
    R_δ    = Parameter(Diagonal(control_params.R_δ  .* dt), m)
    R_Δδ   = Parameter(Diagonal(control_params.R_Δδ ./ dt), m)
    W_β    = Parameter(control_params.W_β .* dt, m)
    W_r    = Parameter(control_params.W_r .* dt, m)
    q_curr = Parameter(Array(qs[1]), m)
    δ_curr = Parameter([us[1].δ], m)
    A      = Parameter{Matrix{T},typeof(identity),true}[]
    B      = Parameter{Matrix{T},typeof(identity),true}[]
    B0     = Parameter{Matrix{T},typeof(identity),true}[]
    Bf     = Parameter{Matrix{T},typeof(identity),true}[]
    c      = Parameter{Vector{T},typeof(identity),true}[]
    H      = Parameter{Matrix{T},typeof(identity),true}[]
    G      = Parameter{Vector{T},typeof(identity),true}[]
    δ_min  = Parameter{Vector{T},typeof(identity),true}[]
    δ_max  = Parameter{Vector{T},typeof(identity),true}[]
    Δδ_min = Parameter{Vector{T},typeof(identity),true}[]
    Δδ_max = Parameter{Vector{T},typeof(identity),true}[]

    q  = [Variable(m) for i in 1:4, t in 1:N_short+N_long+1]
    δ  = [Variable(m) for i in 1:1, t in 1:N_short+N_long+1]
    σ  = [Variable(m) for i in 1:2, t in 1:N_short+N_long]
    Δδ = [Variable(m) for t in 1:N_short+N_long]
    @constraint(m, vec(σ) >= 0)
    @constraint(m, diff(δ[1,:]) == Δδ)

    @constraint(m, q[:,1] == q_curr)
    @constraint(m, δ[:,1] == δ_curr)
    for t in 1:N_short
        ZOHt = linearize(linearize(dynamics, qs[t], [us[t]; ps[t]]),
                         qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1))
        # ZOHt = linearize(dynamics, qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1))
        At = push!(A, Parameter(Array(ZOHt.A), m))[end]
        Bt = push!(B, Parameter(Array(ZOHt.B), m))[end]
        ct = push!(c, Parameter(Array(ZOHt.c), m))[end]
        @constraint(m, At*q[:,t] + Bt*δ[:,t] + ct == q[:,t+1])
    end

    for t in N_short+1:N_short+N_long
        FOHt = linearize(linearize(dynamics, qs[t], [us[t]; ps[t]]),
                         qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1))
        # FOHt = linearize(dynamics, qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1))
        At  = push!(A,  Parameter(Array(FOHt.A),  m))[end]
        B0t = push!(B0, Parameter(Array(FOHt.B0), m))[end]
        Bft = push!(Bf, Parameter(Array(FOHt.Bf), m))[end]
        ct  = push!(c,  Parameter(Array(FOHt.c),  m))[end]
        @constraint(m, At*q[:,t] + B0t*δ[:,t] + Bft*δ[:,t+1] + ct == q[:,t+1])
    end

    δ_hwmax, δ̇_max = dynamics.control_limits.δ_max, control_params.δ̇_max
    for t in 1:N_short+N_long
        Uxt = ps[t+1][1]
        Fxft, Fxrt = longitudinal_tire_forces(dynamics.longitudinal_params, us[t+1][2])
        δ_mint, δ_maxt, Ht, Gt = stable_limits(dynamics.bicycle_model, Uxt, Fxft, Fxrt)
        Ht = push!(H, Parameter(Array(Ht), m))[end]
        Gt = push!(G, Parameter(Array(Gt), m))[end]
        δ_mint = push!(δ_min, Parameter([max(δ_mint, -δ_hwmax)], m))[end]
        δ_maxt = push!(δ_max, Parameter([min(δ_maxt, δ_hwmax)], m))[end]
        Uy_r = q[1:2,t+1]
        σt   = [σ[1,t],σ[1,t],σ[2,t],σ[2,t]]
        Δδt  = [Δδ[t]]
        Δδ_mint = push!(Δδ_min, Parameter([-δ̇_max*dt[t]], m))[end]
        Δδ_maxt = push!(Δδ_max, Parameter([δ̇_max*dt[t]],  m))[end]
        @constraint(m, δ[:,t+1] <= δ_maxt)
        @constraint(m, δ[:,t+1] >= δ_mint)
        @constraint(m, Ht*Uy_r - Gt <= σt)
        @constraint(m, Δδt <= Δδ_maxt)
        @constraint(m, Δδt >= Δδ_mint)
    end

    Δψ   = q[3,2:end]
    e    = q[4,2:end]
    δ⁺   = δ[1,2:end]
    σ1   = σ[1,:]
    σ2   = σ[2,:]
    obj  = @expression transpose(Δψ)*Q_Δψ*Δψ +
                       transpose(e)*Q_e*e +
                       transpose(δ⁺)*R_δ*δ⁺ +
                       transpose(Δδ)*R_Δδ*Δδ +
                       W_β⋅σ1 + W_r⋅σ2
    @objective(m, Minimize, obj)
    m, LateralTrackingQPVariables(q, δ, σ), LateralTrackingQPParams(Q_Δψ, Q_e, R_δ, R_Δδ, W_β, W_r, q_curr, δ_curr,
                                                                    A, B, B0, Bf, c, H, G, δ_min, δ_max, Δδ_min, Δδ_max)
end

function update_QP!(mpc::TrajectoryTrackingMPC, QPP::LateralTrackingQPParams)
    dynamics = mpc.tracking_dynamics
    control_params  = mpc.control_params
    time_steps = mpc.time_steps
    N_short, N_long, dt = time_steps.N_short, time_steps.N_long, time_steps.dt
    qs, us, ps = mpc.qs, mpc.us, mpc.ps

    QPP.Q_Δψ().diag .= control_params.Q_Δψ .* dt
    QPP.Q_e().diag  .= control_params.Q_e  .* dt
    QPP.R_δ().diag  .= control_params.R_δ  .* dt
    QPP.R_Δδ().diag .= control_params.R_Δδ ./ dt
    QPP.W_β()    .= control_params.W_β .* dt
    QPP.W_r()    .= control_params.W_r .* dt
    QPP.q_curr() .= qs[1]
    QPP.δ_curr() .= us[1].δ

    for t in 1:N_short
        ZOHt = linearize(linearize(dynamics, qs[t], [us[t]; ps[t]]),
                         qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1))
        # ZOHt = linearize(dynamics, qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1))
        QPP.A[t]() .= ZOHt.A
        QPP.B[t]() .= ZOHt.B
        QPP.c[t]() .= ZOHt.c
    end
    for t in N_short+1:N_short+N_long
        FOHt = linearize(linearize(dynamics, qs[t], [us[t]; ps[t]]),
                         qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1))
        # FOHt = linearize(dynamics, qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1))
        QPP.A[t]() .= FOHt.A
        QPP.B0[t-N_short]() .= FOHt.B0
        QPP.Bf[t-N_short]() .= FOHt.Bf
        QPP.c[t]() .= FOHt.c
    end
    δ_max, δ̇_max = dynamics.control_limits.δ_max, control_params.δ̇_max
    for t in 1:N_short+N_long
        Uxt = ps[t+1][1]
        Fxft, Fxrt = longitudinal_tire_forces(dynamics.longitudinal_params, us[t+1][2])
        δ_mint, δ_maxt, Ht, Gt = stable_limits(dynamics.bicycle_model, Uxt, Fxft, Fxrt)
        QPP.H[t]() .= Ht
        QPP.G[t]() .= Gt
        QPP.δ_min[t]() .= max(δ_mint, -δ_max)
        QPP.δ_max[t]() .= min(δ_maxt, δ_max)
        QPP.Δδ_min[t]() .= -δ̇_max*dt[t]
        QPP.Δδ_max[t]() .= δ̇_max*dt[t]
    end
end

function get_next_control(mpc::TrajectoryTrackingMPC, variables::LateralTrackingQPVariables)
    BicycleControl(mpc.dynamics.longitudinal_params,
                   BicycleControl2(value(mpc.model, variables.δ[2]), mpc.us[2].Fx))
end
