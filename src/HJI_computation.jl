# export HJICache, GriddedInterpolation, Gridded, Linear, SVector    # TODO: I shouldn't have to do this

@maintain_type struct SimpleCarState{T} <: FieldVector{4,T}
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
    cψ, sψ = sincos(-us.ψ)
    ΔE, ΔN = @SMatrix([cψ sψ; -sψ cψ])*SVector(them.E - us.E, them.N - us.N)
    HJIRelativeState(ΔE, ΔN, adiff(them.ψ, us.ψ), us.Ux, us.Uy, them.V, us.r)
end

struct HJICache
    grid_knots::NTuple{7,Vector{Float32}}
    V ::GriddedInterpolation{Float32,7,Float32,Gridded{Linear},NTuple{7,Vector{Float32}}}
    ∇V::GriddedInterpolation{SVector{7,Float32},7,SVector{7,Float32},Gridded{Linear},NTuple{7,Vector{Float32}}}
end

function placeholder_HJICache()
    grid_knots = tuple((Float32[-1000., 1000.] for i in 1:7)...)
    V  = interpolate(Float32, Float32, grid_knots, zeros(Float32,2,2,2,2,2,2,2), Gridded(Linear()))
    ∇V = interpolate(Float32, SVector{7,Float32}, grid_knots, zeros(SVector{7,Float32},2,2,2,2,2,2,2), Gridded(Linear()))
    HJICache(grid_knots, V, ∇V)
end

function HJICache(fname::String)
    if endswith(fname, ".mat")
        HJIdata = matread(fname)
        grid_knots = tuple((x -> convert(Vector{Float32}, vec(x))).(HJIdata["avoid_set"]["g"]["vs"])...)
        V  = interpolate(Float32, Float32, grid_knots,
                         HJIdata["avoid_set"]["data"][1,1], Gridded(Linear()))
        ∇V = interpolate(Float32, SVector{7,Float32}, grid_knots,
                         SVector.(HJIdata["avoid_set"]["deriv"]...), Gridded(Linear()))
        HJICache(grid_knots, V, ∇V)
    elseif endswith(fname, ".jld2")
        @load fname grid_knots V_raw ∇V_raw
        V  = interpolate(Float32, Float32, grid_knots, V_raw, Gridded(Linear()))
        ∇V = interpolate(Float32, SVector{7,Float32}, grid_knots, Array(reinterpret(SVector{7,Float32}, ∇V_raw)), Gridded(Linear()))
        HJICache(grid_knots, V, ∇V)
    else
        error("Unknown file type for loading HJICache")
    end
end

function save(fname::String, cache::HJICache)
    grid_knots = cache.grid_knots
    V_raw  = cache.V.coefs
    ∇V_raw = Array(reinterpret(Float32, cache.∇V.coefs))
    @save fname grid_knots V_raw ∇V_raw
end

function Base.getindex(cache::HJICache, x::HJIRelativeState{T}) where {T}
    if all(cache.grid_knots[i][1] <= x[i] <= cache.grid_knots[i][end] for i in 1:length(cache.grid_knots))
        (V=cache.V(x[1], x[2], x[3], x[4], x[5], x[6], x[7]), ∇V=cache.∇V(x[1], x[2], x[3], x[4], x[5], x[6], x[7]))    # avoid splatting penalty
    else
        (V=cache.V(max(cache.grid_knots[1][1], min(x[1], cache.grid_knots[1][end])), 
                   max(cache.grid_knots[2][1], min(x[2], cache.grid_knots[2][end])), 
                   max(cache.grid_knots[3][1], min(x[3], cache.grid_knots[3][end])), 
                   max(cache.grid_knots[4][1], min(x[4], cache.grid_knots[4][end])), 
                   max(cache.grid_knots[5][1], min(x[5], cache.grid_knots[5][end])),
                   max(cache.grid_knots[6][1], min(x[6], cache.grid_knots[6][end])),
                   max(cache.grid_knots[7][1], min(x[7], cache.grid_knots[7][end]))), 
        ∇V=cache.∇V(max(cache.grid_knots[1][1], min(x[1], cache.grid_knots[1][end])), 
                    max(cache.grid_knots[2][1], min(x[2], cache.grid_knots[2][end])), 
                    max(cache.grid_knots[3][1], min(x[3], cache.grid_knots[3][end])), 
                    max(cache.grid_knots[4][1], min(x[4], cache.grid_knots[4][end])), 
                    max(cache.grid_knots[5][1], min(x[5], cache.grid_knots[5][end])),
                    max(cache.grid_knots[6][1], min(x[6], cache.grid_knots[6][end])),
                    max(cache.grid_knots[7][1], min(x[7], cache.grid_knots[7][end]))))
    # else
    #     (V=T(Inf), ∇V=zeros(SVector{7,T}))
    end
end

function relative_dynamics(X::VehicleModel, (ΔE, ΔN, Δψ, Ux, Uy, V, r)::StaticVector{7},    # relative state
                                            uR::StaticVector{2},                            # robot control
                                            (ω, a)::StaticVector{2})                        # human control
    bicycle_dynamics = X(BicycleState(ΔE, ΔN, Δψ, Ux, Uy, r), uR)
    sΔψ, cΔψ = sincos(Δψ)
    SVector(
        V*cΔψ - Ux + ΔN*r,
        V*sΔψ - Uy - ΔE*r,
        ω - r,
        bicycle_dynamics[4],
        bicycle_dynamics[5],
        a,
        bicycle_dynamics[6]
    )
end

function optimal_disturbance(X::VehicleModel, relative_state::HJIRelativeState, ∇V::StaticVector{7}, dMode=:min)
    BM, CL = X.bicycle_model, X.control_limits
    μ, m, G = BM.μ, BM.m, BM.G
    κ_max   = CL.κ_max
    Ax_max  = CL.Fx_max/m
    Pmx_max = CL.Px_max/m
    maxA = 0.9*μ*G

    sgn = (dMode == :max ? 1 : -1)

    V = relative_state.V
    lam_w  = ∇V[3]    # psi_rel dot
    lam_Ax = ∇V[6]    # v_dot
    lam_Ay = lam_w/V
    lam_norm = hypot(lam_Ax, lam_Ay)

    if lam_norm < 1e-3
        return SVector{2,Float64}(0, 0)
    else
        desAx = sgn*lam_Ax*maxA/lam_norm
        desAy = sgn*lam_Ay*maxA/lam_norm
        maxAx = min(Ax_max, Pmx_max/V)
        maxAy = κ_max*V*V
        if desAx > maxAx
            if abs(desAy) < maxAy
                maxAy = min(maxAy, sqrt(maxA*maxA - maxAx*maxAx))
            end
            return SVector{2,Float64}(copysign(maxAy, desAy)/V, maxAx)
        else
            if abs(desAy) > maxAy
                if desAx > 0
                    maxAx = min(sqrt(maxA*maxA - maxAy*maxAy), maxAx)
                    return SVector{2,Float64}(copysign(maxAy, desAy)/V, maxAx)
                else
                    return SVector{2,Float64}(copysign(maxAy, desAy)/V, -sqrt(maxA*maxA - maxAy*maxAy))
                end
            else
                return SVector{2,Float64}(desAy/V, maxAx)
            end
        end
    end
end

# for car-car relative dynamics
function optimal_control(X::VehicleModel, relative_state::HJIRelativeState, ∇V::StaticVector{7}, uMode=:max; N=40, M=18)
    BM, LP, CL = X.bicycle_model, X.longitudinal_params, X.control_limits
    L, a, b, h, m, μ, Cαf, Cαr, G, Izz = BM.L, BM.a, BM.b, BM.h, BM.m, BM.μ, BM.Cαf, BM.Cαr, BM.G, BM.Izz
    δ_max, Fx_max, Fx_min = CL.δ_max, CL.Fx_max, CL.Fx_min
    fake_qR = SVector(0., 0., 0., relative_state.Ux, relative_state.Uy, relative_state.r)    # sufficient for lateral_tire_forces

    sgn = (uMode == :max ? 1 : -1)
    V_opt  = -sgn*Inf
    Fx_opt = 0.0
    δ_opt  = 0.0
    for m in 0:(M-1)
        frac_δ = m/(M-1)
        δ = frac_δ * δ_max - (1 - frac_δ) * δ_max
        for n in 0:(N-1)
            frac = n/(N-1)
            Fx = frac*Fx_max + (1 - frac)*Fx_min
            uR = BicycleControl(δ, longitudinal_tire_forces(LP, Fx)...)
            Fyf, Fyr = lateral_tire_forces(BM, fake_qR, uR)
            # Fxf coefficient
            Fxf_coeff = ∇V[4] * cos(δ) / m + ∇V[5] * sin(δ) / m + ∇V[7] * a * sin(δ) / Izz
            # Fxr coefficient
            Fxr_coeff = ∇V[4] / m
            # Fyf coefficient
            Fyf_coeff = -∇V[4] * sin(δ) / m + ∇V[5] * cos(δ) / m + ∇V[7] * a * cos(δ) / Izz
            # Fyr coefficient
            Fyr_coeff = ∇V[5] / m - ∇V[7] * b / Izz 

            V = Fxf_coeff * uR.Fxf + Fyf_coeff * uR.Fxr + Fyf_coeff * Fyf + Fyr_coeff * Fyr

            if (sgn > 0 && V > V_opt) || (sgn < 0 && V < V_opt)
                Fx_opt = Fx
                δ_opt = δ
                V_opt  = V
            end
        end
    end
    BicycleControl2(δ_opt, Fx_opt)
end

# for car-wall relative dynamics
function optimal_control(X::VehicleModel, relative_state::WALLRelativeState, ∇V::StaticVector{5}, uMode=:max; N=40, M=18)
    BM, LP, CL = X.bicycle_model, X.longitudinal_params, X.control_limits
    L, a, b, h, m, μ, Cαf, Cαr, G, Izz = BM.L, BM.a, BM.b, BM.h, BM.m, BM.μ, BM.Cαf, BM.Cαr, BM.G, BM.Izz
    δ_max, Fx_max, Fx_min = CL.δ_max, CL.Fx_max, CL.Fx_min
    fake_qR = SVector(0., 0., 0., relative_state.Ux, relative_state.Uy, relative_state.r)    # sufficient for lateral_tire_forces

    sgn = (uMode == :max ? 1 : -1)
    V_opt  = -sgn*Inf
    Fx_opt = 0.0
    δ_opt  = 0.0
    for m in 0:(M-1)
        frac_δ = m/(M-1)
        δ = frac_δ * δ_max - (1 - frac_δ) * δ_max
        for n in 0:(N-1)
            frac = n/(N-1)
            Fx = frac*Fx_max + (1 - frac)*Fx_min
            uR = BicycleControl(δ, longitudinal_tire_forces(LP, Fx)...)
            Fyf, Fyr = lateral_tire_forces(BM, fake_qR, uR)
            # Fxf coefficient
            Fxf_coeff = ∇V[3] * cos(δ) / m + ∇V[4] * sin(δ) / m + ∇V[5] * a * sin(δ) / Izz
            # Fxr coefficient
            Fxr_coeff = ∇V[3] / m
            # Fyf coefficient
            Fyf_coeff = -∇V[3] * sin(δ) / m + ∇V[4] * cos(δ) / m + ∇V[5] * a * cos(δ) / Izz
            # Fyr coefficient
            Fyr_coeff = ∇V[4] / m - ∇V[5] * b / Izz 

            V = Fxf_coeff * uR.Fxf + Fyf_coeff * uR.Fxr + Fyf_coeff * Fyf + Fyr_coeff * Fyr

            if (sgn > 0 && V > V_opt) || (sgn < 0 && V < V_opt)
                Fx_opt = Fx
                δ_opt = δ
                V_opt  = V
            end
        end
    end
    BicycleControl2(δ_opt, Fx_opt)
end

# for car-car relative dynamics
function compute_reachability_constraint(X::VehicleModel, cache::HJICache, relative_state::HJIRelativeState, ϵ,
                                         uR_lin=optimal_control(X, relative_state, cache[relative_state].∇V))    # definitely not the correct choice...
    V, ∇V = cache[relative_state]
    if V > ϵ 
        (M=SVector{2,Float64}(0, 0), b=1.0)
    else
        uH_opt = optimal_disturbance(X, relative_state, ∇V)
        ∇H_uR = ForwardDiff.gradient(uR -> dot(∇V, relative_dynamics(X, SVector(relative_state), uR, SVector(uH_opt))), SVector(uR_lin))
        (M=∇H_uR, b=dot(∇V, relative_dynamics(X, relative_state, uR_lin, uH_opt)) - dot(∇H_uR, uR_lin)) # so that H = dot(∇V, uR) ≈ ∇H_uR*uR + c
    end
end
