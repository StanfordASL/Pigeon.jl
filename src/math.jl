@inline mod2piF(x::T) where {T} = mod(x, 2*T(pi))
@inline function adiff(x, y)
    d = mod2piF(x - y)
    ifelse(d <= π, d, d - 2*oftype(d, π))
end

cumtrapz(y, x, x0=0) = cumsum(vcat(0, diff(x).*(y[1:end-1] .+ y[2:end])/2)) .+ x0
invcumtrapz(y, x, x0=0) = cumsum(vcat(0, 2*diff(x)./(y[1:end-1] .+ y[2:end]))) .+ x0

@inline function distance2(line::Tuple{<:AbstractVector,<:AbstractVector}, x)
    v = line[2] - line[1]
    λ = clamp(v⋅(x - line[1])/(v⋅v), 0, 1)
    p = (1 - λ)*line[1] + λ*line[2]
    (p - x)⋅(p - x)
end
