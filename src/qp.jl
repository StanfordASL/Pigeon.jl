SimpleQP.Parameter(A::AbstractArray, model) = Parameter(identity, A, model)

struct QP_Parameters{T}
    Q_Δψ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    Q_e   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_δ   ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    R_Δδ  ::Parameter{Diagonal{T,Vector{T}},typeof(identity),true}
    W_β   ::Parameter{Vector{T},typeof(identity),true}
    W_r   ::Parameter{Vector{T},typeof(identity),true}
    curr_q::Parameter{Vector{T},typeof(identity),true}
    curr_δ::Parameter{Vector{T},typeof(identity),true}
    A     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B0    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    B1    ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    c     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    H     ::Vector{Parameter{Matrix{T},typeof(identity),true}}
    G     ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_min ::Vector{Parameter{Vector{T},typeof(identity),true}}
    δ_max ::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_min::Vector{Parameter{Vector{T},typeof(identity),true}}
    Δδ_max::Vector{Parameter{Vector{T},typeof(identity),true}}
end

struct QP_Variables
    q::Matrix{Variable}
    δ::Matrix{Variable}
    σ::Matrix{Variable}
end

function construct_QP(mpc)
    BM = mpc.bicycle_model
    U = mpc.control_params
    N_short, N_long = mpc.N_short, mpc.N_long
    dt = mpc.dt
    qs, us, ps = mpc.qs, mpc.us, mpc.ps
    curr_δ = mpc.current_control.δ

    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.WarmStart(), false)
    m = Model(optimizer)

    Q_Δψ   = Parameter(Diagonal(U.Q_Δψ .* dt), m)
    Q_e    = Parameter(Diagonal(U.Q_e .* dt), m)
    R_δ    = Parameter(Diagonal(U.R_δ .* dt), m)
    R_Δδ   = Parameter(Diagonal(U.R_Δδ ./ dt), m)
    W_β    = Parameter(U.W_β .* dt, m)
    W_r    = Parameter(U.W_r .* dt, m)
    curr_q = Parameter(Array(qs[1]), m)
    curr_δ = Parameter([curr_δ], m)
    A  = Parameter{Matrix{Float64},typeof(identity),true}[]
    B  = Parameter{Matrix{Float64},typeof(identity),true}[]
    B0 = Parameter{Matrix{Float64},typeof(identity),true}[]
    B1 = Parameter{Matrix{Float64},typeof(identity),true}[]
    c  = Parameter{Vector{Float64},typeof(identity),true}[]
    H  = Parameter{Matrix{Float64},typeof(identity),true}[]
    G  = Parameter{Vector{Float64},typeof(identity),true}[]
    δ_min  = Parameter{Vector{Float64},typeof(identity),true}[]
    δ_max  = Parameter{Vector{Float64},typeof(identity),true}[]
    Δδ_min = Parameter{Vector{Float64},typeof(identity),true}[]
    Δδ_max = Parameter{Vector{Float64},typeof(identity),true}[]

    q = [Variable(m) for i in 1:4, t in 1:N_short+N_long+1]
    δ = [Variable(m) for i in 1:1, t in 1:N_short+N_long+1]
    σ = [Variable(m) for i in 1:2, t in 1:N_short+N_long]
    for t in 1:size(σ, 2)
        σt = σ[:,t]
        @constraint(m, σt >= 0.)
    end

    δ⁺ = δ[1,2:end]
    δ⁻ = δ[1,1:end-1]
    Δδ = [Variable(m) for t in 1:N_short+N_long]
    @constraint(m, 1.0*δ⁺ - 1.0*δ⁻ == Δδ)

    q1 = q[:,1]
    δ1 = δ[:,1]
    @constraint(m, q1 == curr_q)
    @constraint(m, δ1 == curr_δ)
    for t in 1:N_short
        At, Bt, ct = ZOH(BM, qs[t], us[t], ps[t], dt[t])
        At = push!(A, Parameter(Array(At), m))[end]
        Bt = push!(B, Parameter(Array(Bt), m))[end]
        ct = push!(c, Parameter(Array(ct), m))[end]
        qt0 = q[:,t]
        qt1 = q[:,t+1]
        δ0  = δ[:,t]
        @constraint(m, At*qt0 + Bt*δ0 + ct == qt1)
    end

    for t in N_short+1:N_short+N_long
        At, B0t, B1t, ct = FOH(BM, qs[t], us[t], ps[t], us[t+1], ps[t+1], dt[t])
        At  = push!(A,  Parameter(Array(At),  m))[end]
        B0t = push!(B0, Parameter(Array(B0t), m))[end]
        B1t = push!(B1, Parameter(Array(B1t), m))[end]
        ct  = push!(c,  Parameter(Array(ct),  m))[end]
        qt0 = q[:,t]
        qt1 = q[:,t+1]
        δ0  = δ[:,t]
        δ1  = δ[:,t+1]
        At, B0t, B1t, ct = A[end], B0[end], B1[end], c[end]
        @constraint(m, At*qt0 + B0t*δ0 + B1t*δ1 + ct == qt1)
    end

    δ_hwmax, δ̇_max = U.δ_hwmax, U.δ̇_max
    for t in 1:N_short+N_long
        δ_mint, δ_maxt, Ht, Gt = stable_limits(BM, qs[t+1], us[t+1], ps[t+1])
        Ht = push!(H, Parameter(Array(Ht), m))[end]
        Gt = push!(G, Parameter(Array(Gt), m))[end]
        δ_mint = push!(δ_min, Parameter([max(δ_mint, -δ_hwmax)], m))[end]
        δ_maxt = push!(δ_max, Parameter([min(δ_maxt, δ_hwmax)], m))[end]
        δ0   = δ[:,t]
        δ1   = δ[:,t+1]
        Uy_r = q[1:2,t+1]
        σt   = [σ[1,t],σ[1,t],σ[2,t],σ[2,t]]
        Δδt  = [Δδ[t]]
        Δδ_mint = push!(Δδ_min, Parameter([-δ̇_max*dt[t]], m))[end]
        Δδ_maxt = push!(Δδ_max, Parameter([δ̇_max*dt[t]],  m))[end]
        @constraint(m, δ1 <= δ_maxt)
        @constraint(m, δ1 >= δ_mint)
        @constraint(m, Ht*Uy_r - Gt <= σt)
        @constraint(m, Δδt <= Δδ_maxt)
        @constraint(m, Δδt >= Δδ_mint)
    end

    Δψ   = q[3,2:end]
    e    = q[4,2:end]
    σ1   = σ[1,:]
    σ2   = σ[2,:]
    obj  = @expression transpose(Δψ)*Q_Δψ*Δψ + transpose(e)*Q_e*e + transpose(δ⁺)*R_δ*δ⁺ + transpose(Δδ)*R_Δδ*Δδ + W_β⋅σ1 + W_r⋅σ2
    @objective(m, Minimize, obj)
    m, QP_Variables(q, δ, σ), QP_Parameters(Q_Δψ, Q_e, R_δ, R_Δδ, W_β, W_r, curr_q, curr_δ, A, B, B0, B1, c, H, G, δ_min, δ_max, Δδ_min, Δδ_max)
end
