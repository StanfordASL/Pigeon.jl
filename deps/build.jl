HJI_cache_filename = joinpath(@__DIR__(), "BicycleCAvoid.jld2")
if !isfile(HJI_cache_filename)
    download("https://github.com/schmrlng/Pigeon.jl/releases/download/v0.2.0/BicycleCAvoid.jld2", HJI_cache_filename)
end