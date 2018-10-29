### /HJI_values
function initialize_HJI_values_marker!(m, cache=X1CMPC.HJI_cache)
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
function update_HJI_values_marker!(m, q, cache=X1CMPC.HJI_cache)
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
function initialize_HJI_contour_marker!(m, cache=X1CMPC.HJI_cache)
    m.header.frame_id = "x1"
    m.type = Marker[:LINE_STRIP]
    m.pose.orientation.w = 1.0
    m.scale.x = 0.1
    m.scale.y = 0.1
    m.scale.z = 0.1
    m.frame_locked = true
    m.color = ColorRGBA(1.0, 1.0, 1.0, 1.0)
end
function update_HJI_contour_marker!(m, q, cache=X1CMPC.HJI_cache)
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