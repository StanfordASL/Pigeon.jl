struct CoupledControlParams{T}
    V_min::T
    V_max::T

    k_V::T    # used only for initial seeding of linearization nodes
    k_s::T

    δ̇_max::T

    Q_Δs::T
    Q_Δψ::T
    Q_e::T
    W_β::T
    W_r::T

    R_δ::T
    R_Δδ::T
    R_Fx::T
    R_ΔFx::T
end
function CoupledControlParams(;V_min=1.0,
                               V_max=15.0,
                               k_V=10/4/100,
                               k_s=10/4/10000,
                               δ̇_max=0.344,
                               Q_Δs=1.0,
                               Q_Δψ=1.0,
                               Q_e=1.0,
                               W_β=50/(10*π/180),
                               W_r=50.0,
                               R_δ=0.0,
                               R_Δδ=0.01,
                               R_Fx=0.0,
                               R_ΔFx=0.5)
    CoupledControlParams(V_min, V_max, k_V, k_s, δ̇_max, Q_Δs, Q_Δψ, Q_e, W_β, W_r, R_δ, R_Δδ, R_Fx, R_ΔFx)
end

function CoupledTrajectoryTrackingMPC(vehicle::Dict{Symbol,T}, trajectory::TrajectoryTube{T}; control_params=CoupledControlParams(),
                                      N_short=10, N_long=20, dt_short=T(0.01), dt_long=T(0.2), use_correction_step=true) where {T}
    dynamics = VehicleModel(vehicle)
    current_state = zeros(BicycleState{T})
    current_control = zeros(BicycleControl{T})
    time_steps = MPCTimeSteps(N_short, N_long, dt_short, dt_long, use_correction_step)

    N = 1 + N_short + N_long
    qs = rand(TrackingBicycleState{T}, N)    # not zeros so that construct_lateral_tracking_QP below works
    us = rand(BicycleControl2{T}, N)
    ps = zeros(TrackingBicycleParams{T}, N)
    tracking_dynamics = VehicleModel(vehicle, TrackingBicycleModel(vehicle))
    model, variables, parameters = construct_coupled_tracking_QP(tracking_dynamics, control_params, time_steps, qs, us, ps)
    TrajectoryTrackingMPC(vehicle, trajectory, dynamics, control_params,
                          current_state, current_control, 0, NaN,
                          time_steps,
                          qs, us, ps,
                          tracking_dynamics, model, variables, parameters, false)
end

function compute_linearization_nodes!(mpc::TrajectoryTrackingMPC{T},
                                      qs::Vector{TrackingBicycleState{T}},
                                      us::Vector{BicycleControl2{T}},
                                      ps::Vector{TrackingBicycleParams{T}}) where {T}
    traj = mpc.trajectory
    X  = mpc.tracking_dynamics
    U  = mpc.control_params
    q0 = mpc.current_state
    u0 = mpc.current_control
    TS = mpc.time_steps
    N_short, N_long = TS.N_short, TS.N_long
    ts, dt, prev_ts = TS.ts, TS.dt, TS.prev_ts

    s0, e0, t0 = path_coordinates(traj, q0)
    tj = traj[s0]
    Δs = s0 - traj(ts[1]).s
    Δψ = adiff(q0.ψ, tj.ψ)
    q  = TrackingBicycleState(Δs, q0.Ux, q0.Uy, q0.r, Δψ, e0)
    u  = BicycleControl2(u0)
    p  = TrackingBicycleParams(tj.V, tj.κ, T(0), T(0))
    if mpc.solved
        qs[1] = q
        us[1] = u
        ps[1] = p
        q_interp = interpolate((prev_ts,), vec(reinterpret(TrackingBicycleState{Float64}, value.(Ref(mpc.model), mpc.variables.q))), Gridded(Linear()))
        u_interp = interpolate((prev_ts,), vec(reinterpret(BicycleControl2{Float64}, value.(Ref(mpc.model), mpc.variables.u))), Gridded(Linear()))
        for i in 2:N_short+N_long+1
            t = ts[i]
            if t < prev_ts[end]
                q = q_interp[t]
                u = u_interp[t] .* mpc.variables.u_normalization
            else
                q = q_interp[prev_ts[end]]
                u = u_interp[prev_ts[end]] .* mpc.variables.u_normalization
            end
            s = traj(t).s + q.Δs
            tj = traj[s]
            p = TrackingBicycleParams(tj.V, tj.κ, T(0), T(0))
            qs[i] = q
            us[i] = u
            ps[i] = p
        end
    else    # steady-state estimates more-or-less the same as in decoupled MPC
        s  = s0
        sΔψ, cΔψ = sincos(Δψ)
        V  = q0.Ux*cΔψ - q0.Uy*sΔψ    # hypot(q0.Ux, q0.Uy)
        β0 = atan(q0.Uy, q0.Ux)
        r0 = q0.r
        δ0 = u0.δ
        Fyf0, Fyr0 = lateral_tire_forces(X.bicycle_model, q0, u0)
        for i in 1:N_short+N_long+1
            τ = (i == N_short+N_long+1 ? dt[i-1] : dt[i])
            tj = traj[s]
            Δs = s - traj(ts[i]).s
            A_des = tj.A + U.k_V*(tj.V - V)/τ + (isnan(mpc.time_offset) ? 0 : -U.k_s*Δs/τ/τ)
            A_des = min(max(A_des, (U.V_min - V)/τ), (U.V_max - V)/τ)
            if i == 1
                q̇ = mpc.dynamics(q0, BicycleControl2(u0), LocalRoadGeometry(tj))
                A = (q̇[4] - q0.r*q0.Uy)*cΔψ - (q̇[5] + q0.r*q0.Ux)*sΔψ
                # A = (q̇[4] - q0.r*q0.Uy)*cos(β0) + (q̇[5] + q0.r*q0.Ux)*sin(β0)
            elseif i <= N_short+1
                est = steady_state_estimates(X, V, A_des, tj.κ, num_iters=1, r=r0, β0=β0, δ0=δ0, Fyf0=Fyf0)
                q = TrackingBicycleState(Δs, q0.Ux, q0.Uy, q0.r, adiff(q0.ψ, tj.ψ), e0)
                u = BicycleControl2(est.δ, est.Fxf + est.Fxr)
                p = TrackingBicycleParams(tj.V, tj.κ, T(0), T(0))
                A = est.A
            else
                est = steady_state_estimates(X, V, A_des, tj.κ)
                q = TrackingBicycleState(Δs, est.Ux, est.Uy, est.r, -est.β, T(0))
                u = BicycleControl2(est.δ, est.Fxf + est.Fxr)
                p = TrackingBicycleParams(tj.V, tj.κ, T(0), T(0))
                A = est.A
            end
            qs[i] = q
            us[i] = u
            ps[i] = p
            i == N_short+N_long+1 && break
            V = V + A*τ
            s = s + V*τ + A*τ*τ/2
        end
    end
end

struct TrackingQPParams{T}
    Q_Δs  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    Q_Δψ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    Q_e   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_δ   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_Δδ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_Fx  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_ΔFx ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    W_β   ::Parameter{Vector{T},typeof(identity),true}
    W_r   ::Parameter{Vector{T},typeof(identity),true}
    q_curr::Parameter{Vector{T},typeof(identity),true}
    u_curr::Parameter{Vector{T},typeof(identity),true}
    A     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B0    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    Bf    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    c     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    H     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    G     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_min ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_max ::Vector{Parameter{Vector{T},typeof(identity),true}}
    Fx_max::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_min::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_max::Vector{Parameter{Vector{T},typeof(identity),true}}
end

struct TrackingQPVariables{T}
    q::Matrix{Variable}
    u::Matrix{Variable}
    σ::Matrix{Variable}
    u_normalization::SVector{2,T}
end

function construct_coupled_tracking_QP(dynamics::VehicleModel{T}, control_params, time_steps, qs, us, ps) where {T}
    N_short, N_long, dt = time_steps.N_short, time_steps.N_long, time_steps.dt
    u_normalization = SVector(dynamics.control_limits.δ_max, max(-dynamics.control_limits.Fx_min, dynamics.control_limits.Fx_max))

    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.WarmStart(), true)
    m = Model{T}(optimizer)

    Q_Δs   = Parameter(Diagonal(control_params.Q_Δs .* dt), m)
    Q_Δψ   = Parameter(Diagonal(control_params.Q_Δψ .* dt), m)
    Q_e    = Parameter(Diagonal(control_params.Q_e  .* dt), m)
    R_δ    = Parameter(Diagonal(control_params.R_δ  .* dt), m)
    R_Δδ   = Parameter(Diagonal(control_params.R_Δδ ./ dt), m)
    R_Fx   = Parameter(Diagonal(control_params.R_Fx  .* dt), m)
    R_ΔFx  = Parameter(Diagonal(control_params.R_ΔFx ./ dt), m)
    W_β    = Parameter(control_params.W_β .* dt, m)
    W_r    = Parameter(control_params.W_r .* dt, m)
    q_curr = Parameter(Array(qs[1]), m)
    u_curr = Parameter(Array(us[1] ./ u_normalization), m)
    A      = Parameter{Matrix{T},typeof(identity),true}[]
    B      = Parameter{Matrix{T},typeof(identity),true}[]
    B0     = Parameter{Matrix{T},typeof(identity),true}[]
    Bf     = Parameter{Matrix{T},typeof(identity),true}[]
    c      = Parameter{Vector{T},typeof(identity),true}[]
    H      = Parameter{Matrix{T},typeof(identity),true}[]
    G      = Parameter{Vector{T},typeof(identity),true}[]
    δ_min  = Parameter{Vector{T},typeof(identity),true}[]
    δ_max  = Parameter{Vector{T},typeof(identity),true}[]
    Fx_max = Parameter{Vector{T},typeof(identity),true}[]
    Δδ_min = Parameter{Vector{T},typeof(identity),true}[]
    Δδ_max = Parameter{Vector{T},typeof(identity),true}[]

    q   = [Variable(m) for i in 1:6, t in 1:N_short+N_long+1]    # (Δs, Ux, Uy, r, Δψ, e)
    u   = [Variable(m) for i in 1:2, t in 1:N_short+N_long+1]    # (δ, Fx)
    σ   = [Variable(m) for i in 1:2, t in 1:N_short+N_long]
    Δδ  = [Variable(m) for t in 1:N_short+N_long]
    ΔFx = [Variable(m) for t in 1:N_short+N_long]
    δ   = u[1,:]
    Fx  = u[2,:]
    Ux  = q[2,:]
    @constraint(m, vec(σ) >= fill(T(0), 2*(N_short+N_long)))
    @constraint(m, diff(δ) == Δδ)
    @constraint(m, diff(Fx) == ΔFx)
    @constraint(m, Ux >= fill(control_params.V_min, N_short+N_long+1))    # may need slacks
    @constraint(m, Ux <= fill(control_params.V_max, N_short+N_long+1))
    @constraint(m, Fx >= fill(dynamics.control_limits.Fx_min / u_normalization[2], N_short+N_long+1))

    @constraint(m, q[:,1] == q_curr)
    @constraint(m, u[:,1] == u_curr)
    for t in 1:N_short
        ZOHt = linearize(dynamics, qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1,2))
        At = push!(A, Parameter(Array(ZOHt.A), m))[end]
        Bt = push!(B, Parameter(Array(ZOHt.B .* u_normalization'), m))[end]
        ct = push!(c, Parameter(Array(ZOHt.c), m))[end]
        @constraint(m, At*q[:,t] + Bt*u[:,t] + ct == q[:,t+1])
    end

    for t in N_short+1:N_short+N_long
        FOHt = linearize(dynamics, qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1,2))
        At  = push!(A,  Parameter(Array(FOHt.A),  m))[end]
        B0t = push!(B0, Parameter(Array(FOHt.B0 .* u_normalization'), m))[end]
        Bft = push!(Bf, Parameter(Array(FOHt.Bf .* u_normalization'), m))[end]
        ct  = push!(c,  Parameter(Array(FOHt.c),  m))[end]
        @constraint(m, At*q[:,t] + B0t*u[:,t] + Bft*u[:,t+1] + ct == q[:,t+1])
    end

    δ_max_, δ̇_max   = dynamics.control_limits.δ_max, control_params.δ̇_max
    Fx_max_, Px_max = dynamics.control_limits.Fx_max, dynamics.control_limits.Px_max
    for t in 1:N_short+N_long
        Uxt = qs[t+1].Ux
        Fxft, Fxrt = longitudinal_tire_forces(dynamics.longitudinal_params, us[t+1].Fx)
        δ_mint, δ_maxt, Ht, Gt = stable_limits(dynamics.bicycle_model, Uxt, Fxft, Fxrt)
        Ht = push!(H, Parameter(Array(Ht), m))[end]
        Gt = push!(G, Parameter(Array(Gt), m))[end]
        δ_mint  = push!(δ_min,  Parameter([max(δ_mint, -δ_max_) / u_normalization[1]], m))[end]
        δ_maxt  = push!(δ_max,  Parameter([min(δ_maxt, δ_max_)  / u_normalization[1]], m))[end]
        Fx_maxt = push!(Fx_max, Parameter([min(Px_max/Uxt, Fx_max_) / u_normalization[2]], m))[end]
        Uy_r = q[3:4,t+1]
        σt   = [σ[1,t],σ[1,t],σ[2,t],σ[2,t]]
        Δδt  = [Δδ[t]]
        Δδ_mint = push!(Δδ_min, Parameter([-δ̇_max*dt[t] / u_normalization[1]], m))[end]
        Δδ_maxt = push!(Δδ_max, Parameter([ δ̇_max*dt[t] / u_normalization[1]], m))[end]
        @constraint(m, [δ[t+1]] <= δ_maxt)
        @constraint(m, [δ[t+1]] >= δ_mint)
        @constraint(m, [Fx[t+1]] <= Fx_maxt)
        @constraint(m, Ht*Uy_r - Gt <= σt)
        @constraint(m, Δδt <= Δδ_maxt)
        @constraint(m, Δδt >= Δδ_mint)
    end

    Δs   = q[1,2:end]
    Δψ   = q[5,2:end]
    e    = q[6,2:end]
    δ⁺   = δ[2:end]
    Fx⁺  = Fx[2:end]
    σ1   = σ[1,:]
    σ2   = σ[2,:]
    obj  = @expression transpose(Δs)*Q_Δs*Δs +
                       transpose(Δψ)*Q_Δψ*Δψ +
                       transpose(e)*Q_e*e +
                       transpose(δ⁺)*R_δ*δ⁺ +
                       transpose(Δδ)*R_Δδ*Δδ +
                       transpose(Fx⁺)*R_Fx*Fx⁺ +
                       transpose(ΔFx)*R_ΔFx*ΔFx +
                       W_β⋅σ1 + W_r⋅σ2
    @objective(m, Minimize, obj)
    m, TrackingQPVariables(q, u, σ, u_normalization), TrackingQPParams(Q_Δs, Q_Δψ, Q_e, R_δ, R_Δδ, R_Fx, R_ΔFx, W_β, W_r,
                                                                       q_curr, u_curr, A, B, B0, Bf, c, H, G,
                                                                       δ_min, δ_max, Fx_max, Δδ_min, Δδ_max)
end

function update_QP!(mpc::TrajectoryTrackingMPC, QPP::TrackingQPParams)
    dynamics = mpc.tracking_dynamics
    control_params  = mpc.control_params
    time_steps = mpc.time_steps
    N_short, N_long, dt = time_steps.N_short, time_steps.N_long, time_steps.dt
    qs, us, ps = mpc.qs, mpc.us, mpc.ps
    u_normalization = mpc.variables.u_normalization

    QPP.Q_Δs().diag  .= control_params.Q_Δs .* dt
    QPP.Q_Δψ().diag  .= control_params.Q_Δψ .* dt
    QPP.Q_e().diag   .= control_params.Q_e  .* dt
    QPP.R_δ().diag   .= control_params.R_δ  .* dt
    QPP.R_Δδ().diag  .= control_params.R_Δδ ./ dt
    QPP.R_Fx().diag  .= control_params.R_Fx  .* dt
    QPP.R_ΔFx().diag .= control_params.R_ΔFx ./ dt
    QPP.W_β()    .= control_params.W_β .* dt
    QPP.W_r()    .= control_params.W_r .* dt
    QPP.q_curr() .= qs[1]
    QPP.u_curr() .= us[1] ./ u_normalization

    for t in 1:N_short
        ZOHt = linearize(dynamics, qs[t], StepControl(dt[t], [us[t]; ps[t]]), keep_control_dims=SVector(1,2))
        QPP.A[t]() .= ZOHt.A
        QPP.B[t]() .= (ZOHt.B .* u_normalization')
        QPP.c[t]() .= ZOHt.c
    end
    for t in N_short+1:N_short+N_long
        FOHt = linearize(dynamics, qs[t], RampControl(dt[t], [us[t]; ps[t]], [us[t+1]; ps[t+1]]), keep_control_dims=SVector(1,2))
        QPP.A[t]() .= FOHt.A
        QPP.B0[t-N_short]() .= (FOHt.B0 .* u_normalization')
        QPP.Bf[t-N_short]() .= (FOHt.Bf .* u_normalization')
        QPP.c[t]() .= FOHt.c
    end
    δ_max, δ̇_max   = dynamics.control_limits.δ_max, control_params.δ̇_max
    Fx_max, Px_max = dynamics.control_limits.Fx_max, dynamics.control_limits.Px_max
    for t in 1:N_short+N_long
        Uxt = qs[t+1].Ux
        Fxft, Fxrt = longitudinal_tire_forces(dynamics.longitudinal_params, us[t+1].Fx)
        δ_mint, δ_maxt, Ht, Gt = stable_limits(dynamics.bicycle_model, Uxt, Fxft, Fxrt)
        QPP.H[t]() .= Ht
        QPP.G[t]() .= Gt
        QPP.δ_min[t]()  .= max(δ_mint, -δ_max) / u_normalization[1]
        QPP.δ_max[t]()  .= min(δ_maxt,  δ_max) / u_normalization[1]
        QPP.Fx_max[t]() .= min(Px_max/Uxt, Fx_max) / u_normalization[2]
        QPP.Δδ_min[t]() .= -δ̇_max*dt[t] / u_normalization[1]
        QPP.Δδ_max[t]() .=  δ̇_max*dt[t] / u_normalization[1]
    end
end

function get_next_control(mpc::TrajectoryTrackingMPC, variables::TrackingQPVariables)
    u_normalization = variables.u_normalization
    δ, Fx = value.(Ref(mpc.model), variables.u[:,2])
    BicycleControl(mpc.dynamics.longitudinal_params,
                   BicycleControl2(δ, Fx) .* u_normalization)
end
