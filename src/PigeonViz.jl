module PigeonViz

using Distributed
using RobotOS
using StaticArrays
using Interpolations
using Contour
using JLD2
using RobotOS
using Rotations

@rosimport std_msgs.msg: ColorRGBA, Header
@rosimport geometry_msgs.msg: Point, TransformStamped
@rosimport visualization_msgs.msg: Marker
@rosimport safe_traffic_weaving.msg: XYThV
@rosimport auto_messages.msg: from_autobox
rostypegen(@__MODULE__)
import .std_msgs.msg: ColorRGBA, Header
import .geometry_msgs.msg: Point, TransformStamped
import .visualization_msgs.msg: Marker
import .safe_traffic_weaving.msg: XYThV
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
struct WALLRelativeState{T} <: FieldVector{5,T}
    ΔN::T
    Δψ::T
    Ux::T
    Uy::T
    r::T
end

function HJIRelativeState(us::BicycleState, them::SimpleCarState)
    cψ, sψ = sincos(-us.ψ)
    ΔE, ΔN = @SMatrix([cψ sψ; -sψ cψ])*SVector(them.E - us.E, them.N - us.N)
    HJIRelativeState(ΔE, ΔN, adiff(them.ψ, us.ψ), us.Ux, us.Uy, them.V, us.r)
end
function WALLRelativeState(us::BicycleState, wall::Vector{Float32})
    a, b, c, ϕ = wall;
    xb = us.E;
    yb = us.N;
    ψ = us.ψ + pi/2;
    y_rel = a * xb + b * yb + c;
    WALLRelativeState(y_rel, adiff(ψ, ϕ), us.Ux, us.Uy, us.r)
end
struct HJICache
    grid_knots::NTuple{7,Vector{Float32}}
    V ::GriddedInterpolation{Float32,7,Float32,Gridded{Linear},NTuple{7,Vector{Float32}}}
    ∇V::GriddedInterpolation{SVector{7,Float32},7,SVector{7,Float32},Gridded{Linear},NTuple{7,Vector{Float32}}}
end
struct WALLCache
    grid_knots::NTuple{5,Vector{Float32}}
    V ::GriddedInterpolation{Float32,5,Float32,Gridded{Linear},NTuple{5,Vector{Float32}}}
    ∇V::GriddedInterpolation{SVector{5,Float32},5,SVector{5,Float32},Gridded{Linear},NTuple{5,Vector{Float32}}}
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
function WALLCache(fname::String)
    if endswith(fname, ".mat")
        HJIdata = matread(fname)
        grid_knots = tuple((x -> convert(Vector{Float32}, vec(x))).(HJIdata["avoid_set"]["g"]["vs"])...)
        V  = interpolate(Float32, Float32, grid_knots,
                         HJIdata["avoid_set"]["data"][1,1], Gridded(Linear()))
        ∇V = interpolate(Float32, SVector{5,Float32}, grid_knots,
                         SVector.(HJIdata["avoid_set"]["deriv"]...), Gridded(Linear()))
        WALLCache(grid_knots, V, ∇V)
    elseif endswith(fname, ".jld2")
        @load fname grid_knots V_raw ∇V_raw
        V  = interpolate(Float32, Float32, grid_knots, V_raw, Gridded(Linear()))
        ∇V = interpolate(Float32, SVector{5,Float32}, grid_knots, Array(reinterpret(SVector{5,Float32}, ∇V_raw)), Gridded(Linear()))
        WALLCache(grid_knots, V, ∇V)
    else
        error("Unknown file type for loading WALLCache")
    end
end
function Base.getindex(cache::HJICache, x::HJIRelativeState{T}) where {T}
    if all(cache.grid_knots[i][1] <= x[i] <= cache.grid_knots[i][end] for i in 1:length(cache.grid_knots))
        (V=cache.V(x[1], x[2], x[3], x[4], x[5], x[6], x[7]), ∇V=cache.∇V(x[1], x[2], x[3], x[4], x[5], x[6], x[7]))    # avoid splatting penalty
    else
        (V=T(Inf), ∇V=zeros(SVector{7,T}))
    end
end
function Base.getindex(cache::WALLCache, x::WALLRelativeState{T}) where {T}
    if all(cache.grid_knots[i][1] <= x[i] <= cache.grid_knots[i][end] for i in 1:length(cache.grid_knots))
        (V=cache.V(x[1], x[2], x[3], x[4], x[5]), ∇V=cache.∇V(x[1], x[2], x[3], x[4], x[5]))    # avoid splatting penalty
    else
        (V=T(Inf), ∇V=zeros(SVector{5,T}))
    end
end


# const cache = HJICache(joinpath(@__DIR__, "../deps/BicycleCAvoid.jld2"))
const cache = HJICache(joinpath(@__DIR__, "13_13_9_9_9_9_9__5.0_CAvoid_Fx_40_delta_18.mat"))
const cache_wall = WALLCache(joinpath(@__DIR__, "21_9_9_5_7__5.0_Wall_Fx_40_delta_18.mat"))


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

function initialize_WALL_values_marker!(m)
    m.header.frame_id = "wall"
    m.type = Marker[:TRIANGLE_LIST]
    m.pose.orientation.w = 1.0
    m.scale.x = 1.0
    m.scale.y = 1.0
    m.scale.z = 1.0
    m.frame_locked = true
    Y = cache_wall.grid_knots[1]
    X = collect(range(-3, stop=5, length=17))
    for i in 1:length(Y)-1, j in 1:length(Y)-1
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

function update_WALL_values_marker!(m, q)
    Y = cache_wall.grid_knots[1]
    X = collect(range(-3, stop=5, length=17))
    ct = 1
    for i in 1:length(X)-1, j in 1:length(Y)-1
        for (x,y) in ((X[i]  , Y[j]  ),
                      (X[i+1], Y[j]  ),
                      (X[i+1], Y[j+1]),
                      (X[i]  , Y[j]  ),
                      (X[i]  , Y[j+1]),
                      (X[i+1], Y[j+1]))
            q = WALLRelativeState{Float32}(y, q[2], q[3], q[4], q[5])
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
const WALL_values_marker = Marker()

initialize_HJI_values_marker!(HJI_values_marker)
initialize_WALL_values_marker!(WALL_values_marker)

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
function initialize_WALL_contour_marker!(m)
    m.header.frame_id = "wall"
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
function update_WALL_contour_marker!(m, q)
    Y = cache.grid_knots[1]
    X = collect(range(-3, stop=5, length=17))
    c = contour(X, Y, [wall_cache[WALLRelativeState{Float32}(y, q[2], q[3], q[4], q[5])].V for x in X, y in Y], Float32(0))
    if isempty(c.lines[1].vertices)
        resize!(m.points, 0)
    else
        m.points = [Point(x, y, -0.05) for (x,y) in c.lines[1].vertices]
    end
end

const HJI_contour_marker = Marker()
const WALL_contour_marker = Marker()

initialize_HJI_contour_marker!(HJI_contour_marker)
initialize_WALL_contour_marker!(WALL_contour_marker)

### /xbox_car/xythv
const other_car_state = fill(zeros(SimpleCarState{Float64}))
function other_car_callback(msg::XYThV)
    other_car_state[] = SimpleCarState(msg.x, msg.y, msg.th - pi/2, msg.v)
end

### /from_autobox
function from_autobox_callback(msg::from_autobox, HJI_values_pub, HJI_contour_pub, WALL_values_pub, WALL_contour_pub, wall)
    try
        current_state = BicycleState(msg.E_m, msg.N_m, msg.psi_rad, msg.ux_mps, msg.uy_mps, msg.r_radps)
        relative_state = HJIRelativeState(current_state, other_car_state[])
        update_HJI_values_marker!(HJI_values_marker, relative_state)
        update_HJI_contour_marker!(HJI_contour_marker, relative_state)


        relative_state = WALLRelativeState(current_state, wall)
        # publish wall frame
        a, b, c, θ = wall
        x_wall = msg.E - a * relative_state.ΔN
        y_wall = msg.N - b * relative_state.ΔN
        world_header = Header()
        world_header.frame_id = "world"
        wall_wrt_world = TransformStamped(world_header, "wall", Transform())
        timestamp = RobotOS.now()
        fill_2D_transform!(wall_wrt_world, x_wall, y_wall, θ, timestamp)
        tf_broadcaster[:sendTransform](wall_wrt_world)

        update_HJI_values_marker!(WALL_values_marker, relative_state)
        update_HJI_contour_marker!(WALL_contour_marker, relative_state)

        publish(WALL_values_pub, WALL_values_marker)
        publish(WALL_contour_pub, WALL_contour_marker)
    catch err
        RobotOS.logwarn("Error: $err\n$(stacktrace(catch_backtrace()))")
    end
end

function fill_2D_transform!(t::TransformStamped, x, y, θ, stamp = RobotOS.now())
    t.header.stamp = stamp
    t.transform.translation.x = x
    t.transform.translation.y = y
    q = Quat(RotZYX(θ, 0, 0))    # tf.transformations[:quaternion_from_euler](0.0, 0.0, θ)
    t.transform.rotation.x = q.x
    t.transform.rotation.y = q.y
    t.transform.rotation.z = q.z
    t.transform.rotation.w = q.w
    t
end

### ROS node init
function start_ROS_node()
    init_node("pigeon_viz", anonymous=false)
    other_car = RobotOS.get_param("human", "/xbox_car")
    roadway = RobotOS.get_param("roadway")
    θ, w = roadway["angle"], roadway["lane_width"]
    x0, y0 = roadway["start_mid"]
    x0 += w * sin(θ)
    y0 -= w * cos(θ)
    a, b = sincos(-θ) # a = -sin(θ), b = cos(θ)
    c = sin(θ) * x0 - y0 * cos(θ)
    wall = SVector{4, Float32}([a, b, c, θ])
    HJI_values_pub = Publisher{Marker}("/HJI_values", queue_size=1)
    HJI_contour_pub = Publisher{Marker}("/HJI_contour", queue_size=1)
    Subscriber{from_autobox}("/from_autobox", from_autobox_callback, (HJI_values_pub, HJI_contour_pub, WALL_values_pub, WALL_contour_pub, wall), queue_size=1)
    Subscriber{XYThV}("$(other_car)/xythv", other_car_callback, queue_size=1)
    @spawn spin()
end

end
