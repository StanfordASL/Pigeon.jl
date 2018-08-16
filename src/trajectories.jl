@maintain_type struct TrajectoryNode{T} <: FieldVector{12,T} t::T; s::T; V::T; A::T; E::T; N::T; ψ::T; κ::T; θ::T; ϕ::T; edge_L::T; edge_R::T end
@maintain_type struct TimeInterpolants{T} <: FieldVector{3,T} s::T; V::T; A::T end
@maintain_type struct SpatialInterpolants{T} <: FieldVector{8,T} E::T; N::T; ψ::T; κ::T; θ::T; ϕ::T; edge_L::T; edge_R::T end

TrajectoryNode(t, ti::TimeInterpolants, si::SpatialInterpolants) = TrajectoryNode(t, ti.s, ti.V, ti.A, si.E, si.N, si.ψ, si.κ, si.θ, si.ϕ, si.edge_L, si.edge_R)
LocalRoadGeometry(tn::TrajectoryNode) = LocalRoadGeometry(tn.ψ, tn.κ, tn.θ, tn.ϕ)

struct TrajectoryTube{T}
    t::Vector{T}         # time (s)
    s::Vector{T}         # longitudinal distance along path (m)
    V::Vector{T}         # longitudinal velocity (m/s)
    A::Vector{T}         # longitudinal acceleration (m/s^2)
    E::Vector{T}         # E coordinate (m)
    N::Vector{T}         # N coordinate (m)
    ψ::Vector{T}         # heading (rad)
    κ::Vector{T}         # curvature (1/m)
    θ::Vector{T}         # grade (rad)
    ϕ::Vector{T}         # bank (rad)
    edge_L::Vector{T}    # left lateral deviation bound (m)
    edge_R::Vector{T}    # right lateral deviation bound (m)

    interp_by_s::GriddedInterpolation{SpatialInterpolants{T}, 1, SpatialInterpolants{T}, Gridded{Linear}, Tuple{Vector{T}}, 0}

    function TrajectoryTube{T}(t::Vector{T}, s::Vector{T}, V::Vector{T}, A::Vector{T},
                               E::Vector{T}, N::Vector{T}, ψ::Vector{T}, κ::Vector{T},
                               θ::Vector{T}, ϕ::Vector{T}, edge_L::Vector{T}, edge_R::Vector{T}) where {T}
        @assert length(t) == length(s) == length(V) == length(A) ==
                length(E) == length(N) == length(ψ) == length(κ) ==
                length(θ) == length(ϕ) == length(edge_L) == length(edge_R)
        interp_by_s = interpolate((s,), vec(reinterpret(SpatialInterpolants{T}, [E N ψ κ θ ϕ edge_L edge_R]')), Gridded(Linear()))
        new(t, s, V, A, E, N, ψ, κ, θ, ϕ, edge_L, edge_R, interp_by_s)
    end
end
function TrajectoryTube(t::Vector{T}, s::Vector{T}, V::Vector{T}, A::Vector{T},
                        E::Vector{T}, N::Vector{T}, ψ::Vector{T}, κ::Vector{T},
                        θ::Vector{T}=zero(t), ϕ::Vector{T}=zero(t),
                        edge_L::Vector{T}=fill(T(4), size(t)), edge_R::Vector{T}=fill(-T(4), size(t))) where {T}
    TrajectoryTube{T}(t, s, V, A, E, N, ψ, κ, θ, ϕ, edge_L, edge_R)
end
Base.length(traj::TrajectoryTube) = length(traj.t)
DifferentialDynamicsModels.duration(traj::TrajectoryTube) = traj.t[end] - traj.t[1]
function (traj::TrajectoryTube)(t)
    i = clamp(searchsortedfirst(traj.t, t, 1, length(traj), Base.Order.ForwardOrdering()) - 1, 1, length(traj)-1)
    A = (traj.V[i+1] - traj.V[i])/(traj.t[i+1] - traj.t[i])    # potentially different from traj.A[i]
    dt = t - traj.t[i]
    ti = TimeInterpolants(traj.s[i] + traj.V[i]*dt + A*dt*dt/2, traj.V[i] + A*dt, A)
    si = traj.interp_by_s[ti.s]
    TrajectoryNode(t, ti, si)
end
function Base.getindex(traj::TrajectoryTube, s)
    i = clamp(searchsortedfirst(traj.s, s, 1, length(traj), Base.Order.ForwardOrdering()) - 1, 1, length(traj)-1)
    A = (traj.V[i+1] - traj.V[i])/(traj.t[i+1] - traj.t[i])    # potentially different from traj.A[i]
    ds = s - traj.s[i]
    if abs(A) < 1e-3 || s > traj.s[end]
        dt = ds/traj.V[i]
    else
        dt = (sqrt(2*A*ds + traj.V[i]^2) - traj.V[i])/A
    end
    t = traj.t[i] + dt
    ti = TimeInterpolants(s, traj.V[i] + A*dt, A)
    si = traj.interp_by_s[s]
    TrajectoryNode(t, ti, si)
end

Base.getindex(traj::TrajectoryTube, i::Integer) = SVector(traj.E[i], traj.N[i])    # this overloading is a bit hacky/brittle
function path_coordinates(traj::TrajectoryTube, x)
    d2min, imin = eltype(x)(Inf), 0
    for i in 1:length(traj)-1
        d2 = distance2((traj[i], traj[i+1]), x)
        if d2 < d2min
            d2min, imin = d2, i
        end
    end
    i = imin
    v = traj[i+1] - traj[i]
    w = x - traj[i]
    ds = sqrt(w⋅w - d2min)
    s = traj.s[i] + ds
    e = sqrt(d2min)*sign(cross(v, w))
    A = (traj.V[i+1] - traj.V[i])/(traj.t[i+1] - traj.t[i])    # potentially different from traj.A[i]
    if abs(A) < 1e-3
        dt = ds/traj.V[i]
    else
        dt = (sqrt(2*A*ds + traj.V[i]^2) - traj.V[i])/A
    end
    t = traj.t[i] + dt
    s, e, t
end
path_coordinates(traj::TrajectoryTube, q::BicycleState) = path_coordinates(traj, SVector(q.E, q.N))

function straight_trajectory(len, vel)
    TrajectoryTube([0., len/vel], # t
                   [0., len],     # s
                   [vel, vel],    # V
                   [0., 0.],      # A
                   [0., 0.],      # E
                   [0., len],     # N
                   [0., 0.],      # ψ
                   [0., 0.])      # κ
end
