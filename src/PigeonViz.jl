module PigeonViz

using Distributed
using RobotOS
using StaticArrays
using Interpolations
using Contour
using JLD2
using RobotOS
@rosimport std_msgs.msg: ColorRGBA
@rosimport geometry_msgs.msg: Point
@rosimport visualization_msgs.msg: Marker
@rosimport asl_prototyping.msg: XYThV
@rosimport auto_messages.msg: from_autobox
rostypegen(@__MODULE__)
import .std_msgs.msg: ColorRGBA
import .geometry_msgs.msg: Point
import .visualization_msgs.msg: Marker
import .asl_prototyping.msg: XYThV
import .auto_messages.msg: from_autobox
import StaticArrays: SUnitRange
import Interpolations: GriddedInterpolation

@inline mod2piF(x::T) where {T} = mod(x, 2*T(pi))
@inline function adiff(x, y)
    d = mod2piF(x - y)
    ifelse(d <= π, d, d - 2*oftype(d, π))
end

### HJICache
struct BicycleState{T} <: FieldVector{6,T}
    E::T    # world frame "x" position of CM
    N::T    # world frame "y" position of CM
    ψ::T    # world frame heading of vehicle
    Ux::T   # body frame longitudinal speed
    Uy::T   # body frame lateral speed
    r::T    # yaw rate (dψ/dt)
end
struct SimpleCarState{T} <: FieldVector{4,T}
    E::T    # world frame "x" position of CM
    N::T    # world frame "y" position of CM
    ψ::T    # world frame heading of vehicle
    V::T    # speed of vehicle (velocity assumed to be in the heading direction)
end
struct HJIRelativeState{T} <: FieldVector{7,T}
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
struct HJICache
    grid_knots::NTuple{7,Vector{Float32}}
    V ::GriddedInterpolation{Float32,7,Float32,Gridded{Linear},NTuple{7,Vector{Float32}},0}
    ∇V::GriddedInterpolation{SVector{7,Float32},7,SVector{7,Float32},Gridded{Linear},NTuple{7,Vector{Float32}},0}
end
function HJICache(fname::String)
    @load fname grid_knots V_raw ∇V_raw
    V  = interpolate(Float32, Float32, grid_knots, V_raw, Gridded(Linear()))
    ∇V = interpolate(Float32, SVector{7,Float32}, grid_knots, Array(reinterpret(SVector{7,Float32}, ∇V_raw)), Gridded(Linear()))
    HJICache(grid_knots, V, ∇V)
end
function Base.getindex(cache::HJICache, x::HJIRelativeState{T}) where {T}
    if all(cache.grid_knots[i][1] <= x[i] <= cache.grid_knots[i][end] for i in 1:length(cache.grid_knots))
        (V=cache.V[x[1], x[2], x[3], x[4], x[5], x[6], x[7]], ∇V=cache.∇V[x[1], x[2], x[3], x[4], x[5], x[6], x[7]])    # avoid splatting penalty
    else
        (V=T(Inf), ∇V=zeros(SVector{7,T}))
    end
end

const cache = HJICache(joinpath(@__DIR__, "BicycleCAvoid.jld2"))

### /HJI_values
function initialize_HJI_values_marker!(m)
    m.header.frame_id = "x1"
    m.type = Marker[:TRIANGLE_LIST]
    m.pose.orientation.w = 1.0
    m.scale.x = 1.0
    m.scale.y = 1.0
    m.scale.z = 1.0
    m.frame_locked = true
    X = cache.grid_knots[1]
    Y = cache.grid_knots[2]
    for i in 1:length(X)-1, j in 1:length(Y)-1
        push!(m.points, Point(X[i]  , Y[j]  , -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
        push!(m.points, Point(X[i+1], Y[j]  , -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
        push!(m.points, Point(X[i+1], Y[j+1], -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
        push!(m.points, Point(X[i]  , Y[j]  , -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
        push!(m.points, Point(X[i]  , Y[j+1], -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
        push!(m.points, Point(X[i+1], Y[j+1], -0.1)); push!(m.colors, ColorRGBA(1.0, 1.0, 1.0, 1.0))
    end
end
function update_HJI_values_marker!(m, q)
    X = cache.grid_knots[1]
    Y = cache.grid_knots[2]
    ct = 1
    for i in 1:length(X)-1, j in 1:length(Y)-1
        for (x,y) in ((X[i]  , Y[j]  ),
                      (X[i+1], Y[j]  ),
                      (X[i+1], Y[j+1]),
                      (X[i]  , Y[j]  ),
                      (X[i]  , Y[j+1]),
                      (X[i+1], Y[j+1]))
            q = HJIRelativeState{Float32}(x, y, q[3], q[4], q[5], q[6], q[7])
            color = m.colors[ct]
            color.r, color.g, color.b = value_to_RGB(cache[q].V)
            ct += 1
        end
    end
end
function value_to_RGB(V, V_lo=-3.0, V_hi=20.0, C_lo = SVector(1.0, 0.5, 0.0), C_hi = SVector(0.0, 0.5, 1.0))
    x = clamp(ifelse(V < 0, 0.5*(V_lo - V)/V_lo, 0.5 + 0.5*V/V_hi), 0.0, 1.0)
    (1 - x)*C_lo + x*C_hi
end

const HJI_values_marker = Marker()
initialize_HJI_values_marker!(HJI_values_marker)

### /HJI_contour
function initialize_HJI_contour_marker!(m)
    m.header.frame_id = "x1"
    m.type = Marker[:LINE_STRIP]
    m.pose.orientation.w = 1.0
    m.scale.x = 0.1
    m.scale.y = 0.1
    m.scale.z = 0.1
    m.frame_locked = true
    m.color = ColorRGBA(1.0, 1.0, 1.0, 1.0)
end
function update_HJI_contour_marker!(m, q)
    X = cache.grid_knots[1]
    Y = cache.grid_knots[2]
    c = contour(X, Y, [cache[HJIRelativeState{Float32}(x, y, q[3], q[4], q[5], q[6], q[7])].V for x in X, y in Y], Float32(0))
    if isempty(c.lines[1].vertices)
        resize!(m.points, 0)
    else
        m.points = [Point(x, y, -0.05) for (x,y) in c.lines[1].vertices]
    end
end

const HJI_contour_marker = Marker()
initialize_HJI_contour_marker!(HJI_contour_marker)

### /xbox_car/xythv
const other_car_state = fill(zeros(SimpleCarState{Float64}))
function other_car_callback(msg::XYThV)
    other_car_state[] = SimpleCarState(msg.x, msg.y, msg.th - pi/2, msg.v)
end

### /from_autobox
function from_autobox_callback(msg::from_autobox, HJI_values_pub, HJI_contour_pub)
    try
        current_state = BicycleState(msg.E_m, msg.N_m, msg.psi_rad, msg.ux_mps, msg.uy_mps, msg.r_radps)
        relative_state = HJIRelativeState(current_state, other_car_state[])
        update_HJI_values_marker!(HJI_values_marker, relative_state)
        update_HJI_contour_marker!(HJI_contour_marker, relative_state)
        publish(HJI_values_pub, HJI_values_marker)
        publish(HJI_contour_pub, HJI_contour_marker)
    catch err
        RobotOS.logwarn("Error: $err\n$(stacktrace(catch_backtrace()))")
    end
end

### ROS node init
function start_ROS_node()
    init_node("pigeon_viz", anonymous=false)
    HJI_values_pub = Publisher{Marker}("/HJI_values", queue_size=1)
    HJI_contour_pub = Publisher{Marker}("/HJI_contour", queue_size=1)
    Subscriber{from_autobox}("/from_autobox", from_autobox_callback, (HJI_values_pub, HJI_contour_pub), queue_size=1)
    Subscriber{XYThV}("/xbox_car/xythv", other_car_callback, queue_size=1)
    @spawn spin()
end

end
