import StaticArrays: SUnitRange

@maintain_type struct BicycleState{T} <: FieldVector{6,T}
    E::T    # world frame "x" position of CM
    N::T    # world frame "y" position of CM
    ψ::T    # world frame heading of vehicle
    Ux::T   # body frame longitudinal speed
    Uy::T   # body frame lateral speed
    r::T    # yaw rate (dψ/dt)
end

@maintain_type struct TrackingBicycleState{T} <: FieldVector{4,T}
    Uy::T   # body frame lateral speed
    r::T    # yaw rate
    Δψ::T   # heading error (w.r.t. nominal trajectory)
    e::T    # lateral error (w.r.t. nominal trajectory)
end

@maintain_type struct TrackingBicycleParams{T} <: FieldVector{4,T}
    Ux::T   # body frame longitudinal speed
    κ::T    # nominal trajectory (local) curvature
    θ::T   # nominal trajectory (local) pitch grade
    ϕ::T   # nominal trajectory (local) roll grade
end

@maintain_type struct BicycleControl{T} <: FieldVector{3,T}
    δ::T    # steering angle
    Fxf::T  # front tire longitudinal force (tire frame)
    Fxr::T  # rear tire longidutinal force (tire frame)
end

struct BicycleModel{T} <: DifferentialDynamics
    # Dimensions
    L::T    # wheelbase
    a::T    # distance from front axle to CG
    b::T    # distance from rear axle to CG
    h::T    # CG height

    # Mass and yaw moment of inertia
    m::T    # vehicle mass
    Izz::T  # yaw moment of inertia

    # Tire model parameters
    μ::T    # coefficient of friction
    Cαf::T  # front tire (pair) cornering stiffness
    Cαr::T  # rear tire (pair) cornering stiffness

    # Longitudinal Drag Force Parameters (FxDrag = Cd0 + Cd1*Ux + Cd2*Ux^2)
    Cd0::T  # rolling resistance
    Cd1::T  # linear drag coefficint
    Cd2::T  # quadratic "aero" drag coefficint
end
BicycleModel(vehicle::Dict) = BicycleModel((vehicle[n] for n in fieldnames(BicycleModel))...)

# Coupled Tire Forces - Simplified Model
@inline function fialatiremodel(α, Cα, μ, Fx, Fz)
    Fmax = μ*Fz
    abs(Fx) >= Fmax ? zero(Fmax) : _fialatiremodel(tan(α), Cα, sqrt(Fmax*Fmax - Fx*Fx))
end

@inline function _fialatiremodel(tanα, Cα, Fymax)
    tanα_slide = 3*Fymax/Cα
    ratio = abs(tanα/tanα_slide)
    if ratio <= 1
        -Cα*tanα*(1 - ratio + ratio*ratio/3)
    else
        -Fymax*sign(tanα)
    end
end

@inline function invfialatiremodel(Fy, Cα, μ, Fx, Fz)
    Fmax = μ*Fz
    Fymax = sqrt(Fmax*Fmax - Fx*Fx)
    atan(_invfialatiremodel(Fy, Cα, Fymax))
end

@inline function _invfialatiremodel(Fy, Cα, Fymax)    # returns tanα
    if abs(Fy) >= Fymax
        -(3*Fymax/Cα)*sign(Fy)
    else
        -(1 + cbrt(abs(Fy)/Fymax - 1))*sign(Fy)
    end
end

@inline function lateral_tire_forces(B::BicycleModel{T}, αf, αr, Fxf, Fxr, sδ, cδ, num_iters=3) where {T}
    L, a, b, h, m, μ, Cαf, Cαr, G = B.L, B.a, B.b, B.h, B.m, B.μ, B.Cαf, B.Cαr, T(GRAVITY)
    Fyf = zero(Fxf)
    Fx = Fxf*cδ - Fyf*sδ + Fxr
    for i in 1:num_iters
        Fzf = (m*G*b - h*Fx)/L
        Fyf = fialatiremodel(αf, Cαf, μ, Fxf, Fzf)
        Fx = Fxf*cδ - Fyf*sδ + Fxr
    end
    Fzr = (m*G*a + h*Fx)/L
    Fyr = fialatiremodel(αr, Cαr, μ, Fxr, Fzr)
    Fyf, Fyr
end

@inline function lateral_tire_forces(B::BicycleModel{T}, q::StaticVector{6}, u::StaticVector{3}, num_iters=3) where {T}
    Ux, Uy, r = q[4], q[5], q[6]
    δ, Fxf, Fxr = u[1], u[2], u[3]
    a, b = B.a, B.b

    sδ, cδ = sincos(δ)
    αf = atan(Uy + a*r, Ux) - δ
    αr = atan(Uy - b*r, Ux)
    lateral_tire_forces(B, αf, αr, Fxf, Fxr, sδ, cδ, num_iters)
end

# Nonlinear Bicycle Model
function (B::BicycleModel{T})(q::StaticVector{6}, u::StaticVector{3}, road::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T}
    E, N, ψ, Ux, Uy, r = q[1], q[2], q[3], q[4], q[5], q[6]
    δ, Fxf, Fxr = u[1], u[2], u[3]
    # ψ, κ, θ, ϕ = road[1], road[2], road[3], road[4]
    L, a, b, h, m, Izz, μ = B.L, B.a, B.b, B.h, B.m, B.Izz, B.μ
    Cαf, Cαr, Cd0, Cd1, Cd2 = B.Cαf, B.Cαr, B.Cd0, B.Cd1, B.Cd2

    sψ, cψ = sincos(ψ)
    sδ, cδ = sincos(δ)
    αf = atan(Uy + a*r, Ux) - δ
    αr = atan(Uy - b*r, Ux)
    Fyf, Fyr = lateral_tire_forces(B, αf, αr, Fxf, Fxr, sδ, cδ)
    Fx_drag = -Cd0 - Ux*(Cd1 + Cd2*Ux)
    Fx_grade = 0    # TODO: figure out how roll/pitch are ordered
    Fy_grade = 0
    F̃xf = Fxf*cδ - Fyf*sδ
    F̃yf = Fyf*cδ + Fxf*sδ
    SVector(
        -Ux*sψ - Uy*cψ, # Ux*cψ - Uy*sψ,  # (ψ measured from N)
        Ux*cψ - Uy*sψ,  # Ux*sψ + Uy*cψ,
        r,
        (F̃xf + Fxr + Fx_drag + Fx_grade)/m + r*Uy,
        (F̃yf + Fyr + Fy_grade)/m - r*Ux,
        (a*F̃yf - b*Fyr)/Izz
    )
end
(B::BicycleModel{T})(q::StaticVector{6}, ur::StaticVector{7}=zeros(LocalRoadGeometry{T})) where {T} = B(q, ur[SUnitRange(1,3)], ur[SUnitRange(4,7)])
(B::BicycleModel{T})(qu::StaticVector{9}, r::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T} = B(qu[SUnitRange(1,6)], qu[SUnitRange(7,9)], r)
(B::BicycleModel)(qur::StaticVector{13}) = B(qur[SUnitRange(1,6)], qur[SUnitRange(7,9)], qur[SUnitRange(10,13)])

# Nonlinear Tracking Bicycle Model
function (B::BicycleModel{T})(q::StaticVector{4}, u::StaticVector{3}, params::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T}
    Uy, r, Δψ, e = q[1], q[2], q[3], q[4]
    δ, Fxf, Fxr = u[1], u[2], u[3]
    Ux, κ, θ, ϕ = params[1], params[2], params[3], params[4]
    L, a, b, h, m, Izz, μ = B.L, B.a, B.b, B.h, B.m, B.Izz, B.μ
    Cαf, Cαr, Cd0, Cd1, Cd2 = B.Cαf, B.Cαr, B.Cd0, B.Cd1, B.Cd2

    L, a, b, h, m, Izz, μ, Cαf, Cαr, G = B.L, B.a, B.b, B.h, B.m, B.Izz, B.μ, B.Cαf, B.Cαr, T(GRAVITY)

    sΔψ, cΔψ = sincos(Δψ)
    sδ, cδ = sincos(δ)
    αf = atan(Uy + a*r, Ux) - δ
    αr = atan(Uy - b*r, Ux)
    Fyf, Fyr = lateral_tire_forces(B, αf, αr, Fxf, Fxr, sδ, cδ)
    Fy_grade = 0
    F̃yf = Fyf*cδ + Fxf*sδ
    SVector(
        (F̃yf + Fyr + Fy_grade)/m - r*Ux,
        (a*F̃yf - b*Fyr)/Izz,
        r - Ux*κ,
        Ux*sΔψ + Uy*cΔψ
    )
end
(B::BicycleModel{T})(q::StaticVector{4}, up::StaticVector{7}=zeros(TrackingBicycleParams{T})) where {T} = B(q, up[SUnitRange(1,3)], up[SUnitRange(4,7)])
(B::BicycleModel{T})(qu::StaticVector{7}, p::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T} = B(qu[SUnitRange(1,4)], qu[SUnitRange(5,7)], p)
(B::BicycleModel)(qup::StaticVector{11}) = B(qup[SUnitRange(1,4)], qup[SUnitRange(5,7)], qup[SUnitRange(8,11)])

# Linearization
function ZOH(B::BicycleModel{T}, q, u, p, dt) where {T}
    qup = [q; u; p]
    J = ForwardDiff.jacobian(B, qup)    # next line could theoretically be folded into this with some extra work/DiffResults magic
    f = B(qup)
    Ainds = @SVector [i for i in 1:4]
    ABinds = @SVector [i for i in 1:5]
    AB = J[:,ABinds]
    c = f - AB*qup[ABinds]
    G = [[AB c]; zeros(SMatrix{2,6,T})]
    expGdt = exp(G*dt)
    expGdt[Ainds, Ainds], expGdt[Ainds, SVector(5)], expGdt[Ainds, 6]
end

function FOH(B::BicycleModel{T}, q, u0, p0, u1, p1, dt) where {T}
    qup = [q; u0; p0]
    J = ForwardDiff.jacobian(B, qup)    # next line could theoretically be folded into this with some extra work/DiffResults magic
    f = B(qup)
    Ainds = @SVector [i for i in 1:4]
    Rinds = @SVector [i for i in 6:11]
    c = f - J*qup
    G = [[J c zeros(SMatrix{4,7,T})]; [zeros(SMatrix{7,12,T}) SMatrix{7,7,T}(I)/dt]; zeros(SMatrix{8,19,T})]    # 1:4 = q, 5:11 = up, 12 = 1, 13:19 = Δup
    expGdt = myexpm(G*dt)
    expGdt[Ainds, Ainds], expGdt[Ainds, SVector(5)] - expGdt[Ainds, SVector(13)], expGdt[Ainds, SVector(13)], expGdt[Ainds, Rinds]*qup[Rinds] + expGdt[Ainds, 8+Rinds]*([u1[SVector(2,3)]; p1] - qup[Rinds]) + expGdt[Ainds, 12]
end

# Ux/r Stability Envelope
@inline function stable_limits(B::BicycleModel{T}, q::StaticVector{4}, u::StaticVector{3}, p::StaticVector{4}) where {T}
    L, a, b, h, m, μ, Cαf, Cαr, G = B.L, B.a, B.b, B.h, B.m, B.μ, B.Cαf, B.Cαr, T(GRAVITY)
    Uy, r, Δψ, e = q[1], q[2], q[3], q[4]
    δ, Fxf, Fxr = u[1], u[2], u[3]
    Ux, κ, θ, ϕ = p[1], p[2], p[3], p[4]

    Fx = Fxf + Fxr
    Fy_grade = 0
    Fzf = (m*G*b - h*Fx)/L
    Fzr = (m*G*a + h*Fx)/L
    Ffmax = μ*Fzf
    Frmax = μ*Fzr
    Fyfmax = abs(Fxf) > Ffmax ? zero(Ffmax) : sqrt(Ffmax*Ffmax - Fxf*Fxf)
    Fyrmax = abs(Fxr) > Frmax ? zero(Frmax) : sqrt(Frmax*Frmax - Fxr*Fxr)
    tanαf_slide = 3*Fyfmax/Cαf
    tanαr_slide = 3*Fyrmax/Cαr
    αf_slide = atan(tanαf_slide)
    αr_slide = atan(tanαr_slide)

    # https://ddl.stanford.edu/sites/default/files/2012_Thesis_Bobier_A_Phase_Portrait_Approach_to_Vehicle_Stabilization_and_Envelope_Control.pdf
    δ_max = atan(L*(μ*G + Fy_grade/m)/(Ux*Ux) - tanαr_slide) + αf_slide    # corresponds to U̇y = 0; αr < 0; αf < 0; δ > 0; r > 0
    δ_min = atan(L*(-μ*G + Fy_grade/m)/(Ux*Ux) + tanαr_slide) - αf_slide   # corresponds to U̇y = 0; αr > 0; αf > 0; δ < 0; r < 0
    rC = (μ*G + Fy_grade/m)/Ux
    UyC = -Ux*tanαr_slide + b*rC
    rD= Ux/L*(tan(αf_slide + δ_max) - tanαr_slide)      # αr > 0; αf > 0; δ > 0;
    UyD = Ux*tanαr_slide + b*rD
    mCD = (rD - rC)/(UyD - UyC)
    rE = Ux/L*(tan(-αf_slide + δ_min) + tanαr_slide)    # αr < 0; αf < 0; δ < 0;
    UyE = -Ux*tanαr_slide + b*rE
    rF = (-μ*G + Fy_grade/m)/Ux
    UyF = Ux*tanαr_slide + b*rF
    mEF = (rF - rE)/(UyF - UyE)

    Hveh = @SMatrix [ 1/Ux -b/Ux;    # β max (≈Uy/Ux)
                     -1/Ux  b/Ux;    # β min (≈Uy/Ux)
                      -mCD     1;    # r max
                       mEF    -1]    # r min
    Gveh = SVector(αr_slide, αr_slide, rC - UyC*mCD, -rF + UyF*mEF)
    δ_min, δ_max, Hveh, Gveh
end

# Vehicle Control Parameters
struct ControlParams{T}
    maxFx::T
    maxPower::T
    fwd_frac::T
    rwd_frac::T
    fwb_frac::T
    rwb_frac::T

    k_V::T
    k_s::T

    δ_hwmax::T
    δ̇_max::T

    Q_Δψ::T
    Q_e::T
    W_β::T
    W_r::T

    R_δ::T
    R_Δδ::T
end
function ControlParams(vehicle::Dict,
                       k_V=GRAVITY/4/100,
                       k_s=GRAVITY/4/10000,
                       δ_hwmax=18*π/180,
                       δ̇_max=0.344,
                       Q_Δψ=1/(10*π/180)^2,
                       Q_e=1.0,
                       W_β=50/(10*π/180),
                       W_r=50.0,
                       R_δ=0.0,
                       R_Δδ=0.01/(10*π/180)^2)
    ControlParams((vehicle[n] for n in fieldnames(ControlParams)[1:6])..., k_V, k_s,
                  δ_hwmax, δ̇_max, Q_Δψ, Q_e, W_β, W_r, R_δ, R_Δδ)
end

# Estimates for Longitudinal Dynamics
function steady_state_estimates(B::BicycleModel{T}, U::ControlParams{T}, κ, A_tan, V, num_iters=4, r=V*κ, β0=T(0), δ0=T(0), Fyf0 = T(0)) where {T}
    # assuming vehicle CM tracks nominal trajectory exactly; constant sideslip w.r.t. nominal trajectory
    L, a, b, h, m, Izz, μ, G = B.L, B.a, B.b, B.h, B.m, B.Izz, B.μ, T(GRAVITY)
    Cαf, Cαr, Cd0, Cd1, Cd2 = B.Cαf, B.Cαr, B.Cd0, B.Cd1, B.Cd2
    maxFx, maxPower, fwd_frac, rwd_frac, fwb_frac, rwb_frac = U.maxFx, U.maxPower, U.fwd_frac, U.rwd_frac, U.fwb_frac, U.rwb_frac

    A_rad = V*V*κ
    A_mag = hypot(A_tan, A_rad)
    A_max = μ*G
    if A_mag > A_max             # nominal trajectory total acceleration exceeds friction limit
        if abs(A_rad) > A_max    # nominal trajectory lateral acceleration exceeds friction limit
            A_rad = A_max*sign(A_rad)
            A_tan = zero(A_rad)
            # error("TODO: case when nominal trajectory lateral acceleration exceeds friction limit")
        else    # preserve radial acceleration for path tracking; reduce intended acceleration along path to compensate
            A_tan = sqrt(A_max*A_max - A_rad*A_rad)*sign(A_tan)
        end
    end
    ṙ = A_tan*κ

    i = 1
    β, δ, Fyf = β0, δ0, Fyf0
    local Ux, Uy, Fxr, Fxf
    while true
        sβ, cβ = sincos(β)
        sδ, cδ = sincos(δ)
        Ux, Uy = V*cβ, V*sβ
        Fx_drag = -Cd0 - Ux*(Cd1 + Cd2*Ux)
        Fx_grade = 0    # TODO: figure out how roll/pitch are ordered
        Fy_grade = 0

        Ax = A_tan*cβ - A_rad*sβ             # Ax = U̇x - r*Uy
        Ay = A_tan*sβ + A_rad*cβ             # Ay = U̇y + r*Ux
        Fx = Ax*m - Fx_drag - Fx_grade       # Fx = F̃xf + Fxr = Fxf*cδ - Fyf*sδ + Fxr
        Fx = min(Fx, min(maxFx, maxPower/Ux)*(rwd_frac + fwd_frac*cδ) - Fyf*sδ)    # braking force saturated by friction, not engine, limits
        Fzr, Fzf = (m*G*a + h*Fx)/L, (m*G*b - h*Fx)/L
        Frmax, Ffmax = μ*Fzr, μ*Fzf

        Fxr = clamp((Fx + Fyf*sδ)*ifelse(Fx > 0, rwd_frac/(rwd_frac + fwd_frac*cδ),
                                                 rwb_frac/(rwb_frac + fwb_frac*cδ)),
                    -Frmax, Frmax)
        Fyrmax = sqrt(Frmax*Frmax - Fxr*Fxr)
        Fyr = (Ay*m - Fy_grade - ṙ*Izz/a)/(1 + b/a)
        Fyr = clamp(Fyr, -Fyrmax, Fyrmax)    # TODO: maybe reduce A_tan in the case that Fyr exceeds Fyrmax
        tanαr = _invfialatiremodel(Fyr, Cαr, Fyrmax)

        F̃xf = clamp(Fx - Fxr, -Ffmax, Ffmax)
        F̃yfmax = sqrt(Ffmax*Ffmax - F̃xf*F̃xf)
        F̃yf = clamp((b*Fyr + ṙ*Izz)/a, -F̃yfmax, F̃yfmax)
        Fxf = F̃xf*cδ + F̃yf*sδ
        Fyf = F̃yf*cδ - F̃xf*sδ
        Fyfmax = sqrt(Ffmax*Ffmax - Fxf*Fxf)
        αf = atan(_invfialatiremodel(Fyf, Cαf, Fyfmax))
        δ = atan(Uy + a*r, Ux) - αf

        if i == num_iters
            Ax = (Fxf*cδ - Fyf*sδ + Fxr + Fx_drag + Fx_grade)/m
            Ay = (Fyf*cδ + Fxf*sδ + Fyr + Fy_grade)/m
            A_tan = Ax*cβ + Ay*sβ
            break
        end
        i = i + 1
        β = atan(tanαr + b*r/Ux)
    end

    Fx_actuated = Fxr + Fxf    # ensure that motor force distribution is exactly respected (e.g., for fwd_frac = 0)
    Fxf = ifelse(Fx_actuated > 0, fwd_frac, fwb_frac)*Fx_actuated
    Fxr = ifelse(Fx_actuated > 0, rwd_frac, rwb_frac)*Fx_actuated
    TrackingBicycleState(Uy, r, -β, 0.0), BicycleControl(δ, Fxf, Fxr), TrackingBicycleParams(Ux, κ, 0.0, 0.0), A_tan
end
