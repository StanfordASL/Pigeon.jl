@maintain_type struct SimpleCarState{T} <: FieldVector{4,T}    # TODO: get from SimpleCarModels.jl
    E::T    # world frame "x" position of CM
    N::T    # world frame "y" position of CM
    ψ::T    # world frame heading of vehicle
    V::T    # speed of vehicle (velocity assumed to be in the heading direction)
end

@maintain_type struct HJIRelativeState{T} <: FieldVector{7,T}
    ΔE::T
    ΔN::T
    Δψ::T
    Ux::T
    Uy::T
    V::T
    r::T
end
function HJIRelativeState(us::BicycleState, them::SimpleCarState)
    cψ, sψ = sincos(us.ψ)
    ΔE, ΔN = @SMatrix([cψ sψ; -sψ cψ])*SVector(them.E - us.E, them.N - us.N)
    HJIRelativeState(ΔE, ΔN, adiff(them.ψ, us.ψ), us.Ux, us.Uy, them.V, us.r)
end

mutable struct MPC{T}
    vehicle::Dict{Symbol,T}
    trajectory::TrajectoryTube{T}
    bicycle_model::BicycleModel{T}
    control_params::ControlParams{T}

    other_car_state::SimpleCarState{T}
    ∇V_cache::Interpolations.GriddedInterpolation{SVector{7,Float64},7,SVector{7,Float32},Gridded{Linear},NTuple{7,Array{Float64,1}},0}
    V_cache::Interpolations.GriddedInterpolation{Float32,7,Float32,Gridded{Linear},NTuple{7,Array{Float64,1}},0}
    reachability_linearization::Tuple{SVector{3,Float64},Float64}
    reachability_ϵ::T

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

        other_car_state = Inf*ones(SimpleCarState{Float64})

        # placeholders to be updated in __init__
        ∇V_cache = interpolate(tuple(([-1000.,1000.] for i in 1:7)...),
                               zeros(SVector{7,Float32},2,2,2,2,2,2,2),
                               Gridded(Linear()))
        V_cache = interpolate(tuple(([-1000.,1000.] for i in 1:7)...),
                              zeros(Float32,2,2,2,2,2,2,2),
                              Gridded(Linear()))
        reachability_linearization = (zeros(SVector{3,T}), T(0))
        reachability_ϵ = T(0.2)

        current_state = zeros(BicycleState{T})
        current_control = zeros(BicycleControl{T})
        N = 1 + N_short + N_long
        ts = T.(1:N) # zeros(T, N) (not zeros so that construct_QP below works)
        dt = diff(ts)
        qs = rand(TrackingBicycleState{T}, N) # not zeros so that construct_QP below works
        us = rand(BicycleControl{T}, N)
        ps = rand(TrackingBicycleParams{T}, N)
        mpc = new(vehicle, trajectory, bicycle_model, control_params,
                  other_car_state, ∇V_cache, V_cache, reachability_linearization, reachability_ϵ,
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
    ts[1:N_short+1] .= t0 .+ dt_short*(0:N_short)
    ts[N_short+2:end] .= t0_long .+ dt_long*(1:N_long)
    for i in 1:N_short+N_long    # long-winded diff!
        dt[i] = ts[i+1] - ts[i]
    end
    ts
end
MPC_time_steps!(mpc::MPC, t0) = MPC_time_steps!(mpc.ts, mpc.dt, t0, mpc.N_short, mpc.N_long, mpc.dt_short, mpc.dt_long, mpc.use_correction_step)
MPC_time_steps(t0, N_short=10, N_long=20, dt_short=0.01, dt_long=0.2, use_correction_step=true) =
    MPC_time_steps!(zeros(typeof(t0), 1+N_short+N_long), zeros(typeof(t0), N_short+N_long), t0, N_short, N_long, dt_short, dt_long, use_correction_step)


function in_cache(cache, x)
    all(cache.knots[i][1] < x[i] < cache.knots[i][end] for i in 1:length(cache.knots))
end

function compute_reachability_constraint!(mpc)
    relative_state = HJIRelativeState(mpc.current_state, mpc.other_car_state)
    !in_cache(mpc.∇V_cache, relative_state) && return
    ∇V  = mpc.∇V_cache[relative_state...]
    uH  = optimal_disturbance(mpc.vehicle, relative_state, ∇V)
    uR0 = mpc.current_control
    ∇H_uR = ForwardDiff.gradient(uR -> dot(∇V, relative_dynamics(mpc.bicycle_model,
                                                                 SVector(relative_state),
                                                                 SVector(uR),
                                                                 SVector(uH))), SVector(uR0))
    c = dot(∇V, relative_dynamics(mpc.bicycle_model, relative_state, uR0, uH)) - dot(∇H_uR, uR0)
    # so that dot(uR) ≈ ∇H_uR*uR + c
    mpc.reachability_linearization = (∇H_uR, c - 0.25)
end

function compute_linearization_nodes!(mpc::MPC)
    traj = mpc.trajectory
    B = mpc.bicycle_model
    U = mpc.control_params
    q0 = mpc.current_state
    u0 = mpc.current_control
    N_short, N_long = mpc.N_short, mpc.N_long
    ts, dt = mpc.ts, mpc.dt

    s, e0, t0 = path_coordinates(traj, SVector(q0.E, q0.N))

    V = hypot(q0.Ux, q0.Uy)
    β0 = atan(q0.Uy, q0.Ux)
    r0 = q0.r
    δ0 = u0.δ
    Fyf0, Fyr0 = lateral_tire_forces(B, q0, u0)
    A_des_min, A_des_max = -Inf, Inf
    relative_state = HJIRelativeState(mpc.current_state, mpc.other_car_state)
    apply_reachability_constraint = in_cache(mpc.V_cache, relative_state) && (mpc.V_cache[relative_state...] < mpc.reachability_ϵ)
    M, b = mpc.reachability_linearization
    for i in 1:N_short+N_long+1
        τ = (i == N_short+N_long+1 ? dt[i-1] : dt[i])
        tj = traj[s]
        κ = tj.κ
        A_des = tj.A + U.k_V*(tj.V - V)/τ + (isnan(mpc.time_offset) ? 0 : U.k_s*(traj(ts[i]).s - s)/τ/τ)   # alternatively, - U.k_t*(tj.t - ts[i])/τ/τ
        A_des = min(max(A_des, (U.V_min - V)/τ), (U.V_max - V)/τ)
        if i <= N_short+1
            q = TrackingBicycleState(q0.Uy, q0.r, adiff(q0.ψ, tj.ψ), e0)    # to match paper should be (..., 0.0, 0.0)
            _, u, p, A = steady_state_estimates(B, U, κ, A_des, V, 1, r0, β0, δ0, Fyf0)
            A_des, δ_des = A, u.δ
            if apply_reachability_constraint
                Fx_drag = B.Cd0 + q0.Ux*(B.Cd1 + B.Cd2*q0.Ux)
                Fx_des = A_des*B.m + Fx_drag
                if Fx_des > 0
                    Fxf_des, Fxr_des = U.fwd_frac*Fx_des, U.rwd_frac*Fx_des
                else
                    Fxf_des, Fxr_des = U.fwb_frac*Fx_des, U.rwb_frac*Fx_des
                end
                Hderiv = dot(BicycleControl(δ_des, Fxf_des, Fxr_des), M) + b
                if Hderiv < 0
                    M_Fx = Fx_des < 0 ? U.fwb_frac*M[2] + U.rwb_frac*M[3] :
                                        U.fwd_frac*M[2] + U.rwd_frac*M[3]
                    ΔFx_frac = 5000*abs(M_Fx)        # super hack
                    Δδ_frac  = 10*pi/180*abs(M[1])   # super hack
                    frac_sum = ΔFx_frac + Δδ_frac
                    ΔFx_frac, Δδ_frac = ΔFx_frac/frac_sum, Δδ_frac/frac_sum

                    Fx_des = Fx_des + ΔFx_frac*(-Hderiv)/M_Fx
                    δ_des  =  δ_des +  Δδ_frac*(-Hderiv)/M[1]
                    Fx_des = clamp(Fx_des, mpc.vehicle[:minFx], mpc.vehicle[:maxFx])
                    δ_des  = clamp(δ_des, -mpc.vehicle[:d_max], mpc.vehicle[:d_max])
                    A = (Fx_des - Fx_drag)/B.m
                    Fxf, Fxr = Fx_des < 0 ? (U.fwb_frac*Fx_des, U.rwb_frac*Fx_des) :
                                            (U.fwd_frac*Fx_des, U.rwd_frac*Fx_des)
                    u = BicycleControl(δ_des, Fxf, Fxr)
                end
            end
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

    relative_state = HJIRelativeState(mpc.current_state, mpc.other_car_state)
    apply_reachability_constraint = in_cache(mpc.V_cache, relative_state) && (mpc.V_cache[relative_state...] < mpc.reachability_ϵ)
    for t in 1:N_short
        At, Bt, ct = ZOH(BM, qs[t], us[t], ps[t], dt[t])
        QPP.A[t]() .= At
        QPP.B[t]() .= Bt
        QPP.c[t]() .= ct
        if apply_reachability_constraint
            M, b = mpc.reachability_linearization
            QPP.Mδ_HJI[t]() .= M[1]
            QPP.bδ_HJI[t]() .= b + M[2]*us[t+1][2] + M[3]*us[t+1][3]
        else
            QPP.Mδ_HJI[t]() .= 0
            QPP.bδ_HJI[t]() .= 1
        end
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
    xs = []
    us = []
    ps = []
    for t in 0:dt:mpc.trajectory.t[end]
        push!(qs, mpc.current_state)
        push!(us, mpc.current_control)
        MPC_time_steps!(mpc, t)
        compute_linearization_nodes!(mpc)
        update_QP!(mpc)
        solve!(mpc.model)
        mpc.current_state = propagate(mpc.bicycle_model, mpc.current_state, StepControl(dt, mpc.current_control))
        mpc.current_control = BicycleControl(value(mpc.model, mpc.variables.δ[2]), mpc.us[2][2], mpc.us[2][3])
        push!(xs, mpc.qs[1])
        push!(ps, mpc.ps[1])
    end
    qs, xs, us, ps
end
