import StaticArrays: SUnitRange

@maintain_type struct BicycleState{T} <: FieldVector{6,T}
    E::T    # world frame "x" position of CM
    N::T    # world frame "y" position of CM
    ψ::T    # world frame heading of vehicle
    Ux::T   # body frame longitudinal speed
    Uy::T   # body frame lateral speed
    r::T    # yaw rate (dψ/dt)
end

@maintain_type struct LocalRoadGeometry{T} <: FieldVector{4,T}
    ψ::T    # nominal trajectory heading
    κ::T    # nominal trajectory (local) curvature
    θ::T    # nominal trajectory (local) pitch grade
    ϕ::T    # nominal trajectory (local) roll grade
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
    θ::T    # nominal trajectory (local) pitch grade
    ϕ::T    # nominal trajectory (local) roll grade
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
    G::T    # gravitational acceleration
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
    F_max = μ*Fz
    abs(Fx) >= F_max ? zero(F_max) : _fialatiremodel(tan(α), Cα, sqrt(F_max*F_max - Fx*Fx))
end

@inline function _fialatiremodel(tanα, Cα, Fy_max)
    tanα_slide = 3*Fy_max/Cα
    ratio = abs(tanα/tanα_slide)
    if ratio <= 1
        -Cα*tanα*(1 - ratio + ratio*ratio/3)
    else
        -Fy_max*sign(tanα)
    end
end

@inline function invfialatiremodel(Fy, Cα, μ, Fx, Fz)
    F_max = μ*Fz
    Fy_max = sqrt(F_max*F_max - Fx*Fx)
    atan(_invfialatiremodel(Fy, Cα, Fy_max))
end

@inline function _invfialatiremodel(Fy, Cα, Fy_max)    # returns tanα
    if abs(Fy) >= Fy_max
        -(3*Fy_max/Cα)*sign(Fy)
    else
        -(1 + cbrt(abs(Fy)/Fy_max - 1))*sign(Fy)
    end
end

@inline function lateral_tire_forces(B::BicycleModel, αf, αr, Fxf, Fxr, sδ, cδ, num_iters=3)
    L, a, b, h, m, μ, Cαf, Cαr, G = B.L, B.a, B.b, B.h, B.m, B.μ, B.Cαf, B.Cαr, B.G
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

@inline function lateral_tire_forces(B::BicycleModel, q::StaticVector{6}, u::StaticVector{3}, num_iters=3)
    Ux, Uy, r = q[4], q[5], q[6]
    δ, Fxf, Fxr = u[1], u[2], u[3]
    a, b = B.a, B.b

    sδ, cδ = sincos(δ)
    αf = atan(Uy + a*r, Ux) - δ
    αr = atan(Uy - b*r, Ux)
    lateral_tire_forces(B, αf, αr, Fxf, Fxr, sδ, cδ, num_iters)
end

# Bicycle Model
function (B::BicycleModel{T})((E, N, ψ, Ux, Uy, r)::StaticVector{6},
                              (δ, Fxf, Fxr)::StaticVector{3},
                              (ψᵣ, κ, θ, ϕ)::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T}
    a, b, m, Izz, Cd0, Cd1, Cd2 = B.a, B.b, B.m, B.Izz, B.Cd0, B.Cd1, B.Cd2

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
(B::BicycleModel{T})(q::StaticVector{6}, ur::StaticVector{7}) where {T} = B(q, ur[SUnitRange(1,3)], ur[SUnitRange(4,7)])
(B::BicycleModel{T})(qu::StaticVector{9}, r::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T} = B(qu[SUnitRange(1,6)], qu[SUnitRange(7,9)], r)
(B::BicycleModel)(qur::StaticVector{13}) = B(qur[SUnitRange(1,6)], qur[SUnitRange(7,9)], qur[SUnitRange(10,13)])

# Lateral Tracking Bicycle Model
function (B::BicycleModel{T})((Uy, r, Δψ, e)::StaticVector{4},
                              (δ, Fxf, Fxr)::StaticVector{3},
                              (Ux, κ, θ, ϕ)::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T}
    a, b, m, Izz = B.a, B.b, B.m, B.Izz

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
(B::BicycleModel{T})(q::StaticVector{4}, up::StaticVector{7}) where {T} = B(q, up[SUnitRange(1,3)], up[SUnitRange(4,7)])
(B::BicycleModel{T})(qu::StaticVector{7}, p::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T} = B(qu[SUnitRange(1,4)], qu[SUnitRange(5,7)], p)
(B::BicycleModel)(qup::StaticVector{11}) = B(qup[SUnitRange(1,4)], qup[SUnitRange(5,7)], qup[SUnitRange(8,11)])

# Uy/r Stability Envelope
@inline function stable_limits(B::BicycleModel, Ux, Fxf, Fxr)
    L, a, b, h, m, μ, Cαf, Cαr, G = B.L, B.a, B.b, B.h, B.m, B.μ, B.Cαf, B.Cαr, B.G

    Fx = Fxf + Fxr
    Fy_grade = 0
    Fzf = (m*G*b - h*Fx)/L
    Fzr = (m*G*a + h*Fx)/L
    Ff_max = μ*Fzf
    Fr_max = μ*Fzr
    Fyf_max = abs(Fxf) > Ff_max ? zero(Ff_max) : sqrt(Ff_max*Ff_max - Fxf*Fxf)
    Fyr_max = abs(Fxr) > Fr_max ? zero(Fr_max) : sqrt(Fr_max*Fr_max - Fxr*Fxr)
    tanαf_slide = 3*Fyf_max/Cαf
    tanαr_slide = 3*Fyr_max/Cαr
    αf_slide = atan(tanαf_slide)
    αr_slide = atan(tanαr_slide)

    # https://ddl.stanford.edu/sites/default/files/publications/2012_Thesis_Bobier_A_Phase_Portrait_Approach_to_Vehicle_Stabilization_and_Envelope_Control.pdf
    δ_max = atan(L*( μ*G + Fy_grade/m)/(Ux*Ux) - tanαr_slide) + αf_slide    # corresponds to U̇y = 0; αr < 0; αf < 0; δ > 0; r > 0
    δ_min = atan(L*(-μ*G + Fy_grade/m)/(Ux*Ux) + tanαr_slide) - αf_slide    # corresponds to U̇y = 0; αr > 0; αf > 0; δ < 0; r < 0
    rC  = ( μ*G + Fy_grade/m)/Ux
    UyC = -Ux*tanαr_slide + b*rC
    rD  = Ux/L*(tan( αf_slide + δ_max) - tanαr_slide)     # αr > 0; αf > 0; δ > 0;
    UyD = Ux*tanαr_slide + b*rD
    mCD = (rD - rC)/(UyD - UyC)
    rE  = Ux/L*(tan(-αf_slide + δ_min) + tanαr_slide)    # αr < 0; αf < 0; δ < 0;
    UyE = -Ux*tanαr_slide + b*rE
    rF  = (-μ*G + Fy_grade/m)/Ux
    UyF = Ux*tanαr_slide + b*rF
    mEF = (rF - rE)/(UyF - UyE)

    H_veh = @SMatrix [ 1/Ux -b/Ux;    # β max (≈Uy/Ux)
                      -1/Ux  b/Ux;    # β min (≈Uy/Ux)
                       -mCD     1;    # r max
                        mEF    -1]    # r min
    G_veh = SVector(αr_slide, αr_slide, rC - UyC*mCD, -rF + UyF*mEF)
    δ_min, δ_max, H_veh, G_veh
end

# Bicycle Model with Combined Longitudinal Force Input (Fx = Fxf + Fxr)
@maintain_type struct BicycleControl2{T} <: FieldVector{2,T}
    δ::T    # steering angle
    Fx::T   # combined front/rear tire longitudinal force (sum of each in respective tire frames; approximately, total is applied longitudinally in body frame)
end
BicycleControl2((δ, Fxf, Fxr)::StaticVector{3}) = BicycleControl2(δ, Fxf + Fxr)

struct LongitudinalActuationParams{T}
    fwd_frac::T
    rwd_frac::T
    fwb_frac::T
    rwb_frac::T
end
LongitudinalActuationParams(vehicle::Dict) = LongitudinalActuationParams((vehicle[n] for n in fieldnames(LongitudinalActuationParams))...)
function longitudinal_tire_forces(LP::LongitudinalActuationParams, Fx)
    Fx > 0 ? (Fx*LP.fwd_frac, Fx*LP.rwd_frac) :
             (Fx*LP.fwb_frac, Fx*LP.rwb_frac)
end
BicycleControl(LP::LongitudinalActuationParams, (δ, Fx)::StaticVector{2}) = BicycleControl(δ, longitudinal_tire_forces(LP, Fx)...)

struct ControlLimits{T}
    Fx_max::T
    Fx_min::T
    Px_max::T
    δ_max::T
    κ_max::T
end
ControlLimits(vehicle::Dict) = ControlLimits((vehicle[n] for n in fieldnames(ControlLimits))...)
function apply_control_limits(CL::ControlLimits, (δ, Fx)::StaticVector{2}, Ux)
    Fx_max, Fx_min, Px_max, δ_max = CL.Fx_max, CL.Fx_min, CL.Px_max, CL.δ_max
    BicycleControl2(clamp(δ, -δ_max, δ_max),
                    max(min(Fx, Fx_max, Px_max/Ux), Fx_min))
end

struct VehicleModel{T} <: DifferentialDynamics
    bicycle_model::BicycleModel{T}
    longitudinal_params::LongitudinalActuationParams{T}
    control_limits::ControlLimits{T}
end
VehicleModel(vehicle::Dict) = VehicleModel(BicycleModel(vehicle), LongitudinalActuationParams(vehicle), ControlLimits(vehicle))

## Full Bicycle Dynamics
(X::VehicleModel{T})(q::StaticVector{6}, u::StaticVector{2}, p::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T} =
    X.bicycle_model(q, BicycleControl(X.longitudinal_params, apply_control_limits(X.control_limits, u, q[4])), p)
(X::VehicleModel{T})(q::StaticVector{6}, ur::StaticVector{6}) where {T} = X(q, ur[SUnitRange(1,2)], ur[SUnitRange(3,6)])
(X::VehicleModel{T})(qu::StaticVector{8}, r::StaticVector{4}=zeros(LocalRoadGeometry{T})) where {T} = X(qu[SUnitRange(1,6)], qu[SUnitRange(7,8)], r)
(X::VehicleModel)(qur::StaticVector{12}) = X(qur[SUnitRange(1,6)], qur[SUnitRange(7,8)], qur[SUnitRange(9,12)])
## Tracking Bicycle Dynamics
(X::VehicleModel{T})(q::StaticVector{4}, u::StaticVector{2}, p::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T} =
    X.bicycle_model(q, BicycleControl(X.longitudinal_params, apply_control_limits(X.control_limits, u, p[1])), p)
(X::VehicleModel{T})(q::StaticVector{4}, up::StaticVector{6}) where {T} = X(q, up[SUnitRange(1,2)], up[SUnitRange(3,6)])
(X::VehicleModel{T})(qu::StaticVector{6}, p::StaticVector{4}=zeros(TrackingBicycleParams{T})) where {T} = X(qu[SUnitRange(1,4)], qu[SUnitRange(5,6)], p)
(X::VehicleModel)(qup::StaticVector{10}) = X(qup[SUnitRange(1,4)], qup[SUnitRange(5,6)], qup[SUnitRange(7,10)])

# Estimates for Longitudinal Dynamics
function steady_state_estimates(X::VehicleModel{T}, V, A_tan, κ;
                                num_iters=4, r=V*κ, β0=T(0), δ0=T(0), Fyf0=T(0)) where {T}
    # assuming vehicle CM tracks nominal trajectory exactly; constant sideslip w.r.t. nominal trajectory
    B, LP, CL = X.bicycle_model, X.longitudinal_params, X.control_limits
    L, a, b, h, m, Izz, μ, G = B.L, B.a, B.b, B.h, B.m, B.Izz, B.μ, B.G
    Cαf, Cαr, Cd0, Cd1, Cd2 = B.Cαf, B.Cαr, B.Cd0, B.Cd1, B.Cd2
    fwd_frac, rwd_frac, fwb_frac, rwb_frac = LP.fwd_frac, LP.rwd_frac, LP.fwb_frac, LP.rwb_frac
    Fx_max, Px_max = CL.Fx_max, CL.Px_max

    A_rad = V*V*κ
    A_mag = hypot(A_tan, A_rad)
    A_max = μ*G
    if A_mag > A_max             # nominal trajectory total acceleration exceeds friction limit
        if abs(A_rad) > A_max    # nominal trajectory lateral acceleration exceeds friction limit
            A_rad = A_max*sign(A_rad)
            A_tan = zero(A_rad)
            # error("TODO: case when nominal trajectory lateral acceleration exceeds friction limit")
        else    # prioritize radial acceleration for path tracking; reduce intended acceleration along path to compensate
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
        Fx = min(Fx, min(Fx_max, Px_max/Ux)*(rwd_frac + fwd_frac*cδ) - Fyf*sδ)    # braking force saturated by friction, not engine, limits
        Fzr, Fzf = (m*G*a + h*Fx)/L, (m*G*b - h*Fx)/L
        Fr_max, Ff_max = μ*Fzr, μ*Fzf

        Fxr = clamp((Fx + Fyf*sδ)*ifelse(Fx > 0, rwd_frac/(rwd_frac + fwd_frac*cδ),
                                                 rwb_frac/(rwb_frac + fwb_frac*cδ)),
                    -Fr_max, Fr_max)
        Fyr_max = sqrt(Fr_max*Fr_max - Fxr*Fxr)
        Fyr = (Ay*m - Fy_grade - ṙ*Izz/a)/(1 + b/a)
        Fyr = clamp(Fyr, -Fyr_max, Fyr_max)    # TODO: maybe reduce A_tan in the case that Fyr exceeds Fyr_max
        tanαr = _invfialatiremodel(Fyr, Cαr, Fyr_max)

        F̃xf = clamp(Fx - Fxr, -Ff_max, Ff_max)
        F̃yfmax = sqrt(Ff_max*Ff_max - F̃xf*F̃xf)
        F̃yf = clamp((b*Fyr + ṙ*Izz)/a, -F̃yfmax, F̃yfmax)
        Fxf = F̃xf*cδ + F̃yf*sδ
        Fyf = F̃yf*cδ - F̃xf*sδ
        Fyf_max = sqrt(Ff_max*Ff_max - Fxf*Fxf)
        αf = atan(_invfialatiremodel(Fyf, Cαf, Fyf_max))
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

    TrackingBicycleState(Uy, r, -β, 0.0), BicycleControl2(δ, Fxf + Fxr), TrackingBicycleParams(Ux, κ, 0.0, 0.0), A_tan
end
