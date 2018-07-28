const GRAVITY = 9.80665

@inline mod2piF(x::T) where {T<:AbstractFloat} = mod(x, 2*T(pi))
@inline function adiff(x::T, y::T) where {T<:AbstractFloat}
    d = mod2piF(x - y)
    d <= π ? d : d - 2*T(π)
end

cumtrapz(y::Vector, x::Vector, x0=0) = cumsum(vcat(0, diff(x).*(y[1:end-1] .+ y[2:end])/2)) .+ x0
invcumtrapz(y::Vector, x::Vector, x0=0) = cumsum(vcat(0, 2*diff(x)./(y[1:end-1] .+ y[2:end]))) .+ x0

@inline function distance2(line::Tuple{<:AbstractVector,<:AbstractVector}, x)
    v = line[2] - line[1]
    λ = clamp(v'*(x - line[1])/(v'*v), 0, 1)
    p = (1 - λ)*line[1] + λ*line[2]
    (p - x)'*(p - x)
end

macro evalmatrixpoly!(B, C, d, Apow_p...)    # probably no reason for this to be a macro
    Apow = Apow_p[1:d]
    p = Apow_p[d+1:end]
    partitioned_inds = reverse(collect(Iterators.partition(1:length(p), d)))
    ex = :(B = fill!($(esc(B)), 0))
    gemm_ex = :(BLAS.gemm!('N','N',one($(esc(p[1]))),$(esc(Apow[d])),B,one($(esc(p[1]))),C))
    highest_order_terms = true
    for inds in partitioned_inds
        ex = quote
            $ex
            C = fill!($(esc(C)), 0)
            for k in diagind(C)
                C[k] = $(esc(p[inds[1]]))
            end
            $(begin
                subex = :(C)
                for i in 2:length(inds)
                    subex = :(BLAS.axpy!($(esc(p[inds[i]])), $(esc(Apow[i-1])), $subex))
                end
                subex
            end)
            $(if !highest_order_terms
                gemm_ex
            end)
            copyto!(B, C)
        end
        highest_order_terms = false
    end
    ex
end

import LinearAlgebra: BlasInt, BlasFloat, checksquare, SingularException

function mygesv!(A::StridedMatrix{T}, B::StridedMatrix{T}, ipiv::Vector{Int}, perm::Vector{Int}) where T
    # LU factorization
    n = checksquare(A)
    @inbounds begin
        for k = 1:n
            # find index max
            kp = k
            amax = abs(zero(T))
            for i = k:n
                absi = abs(A[i,k])
                if absi > amax
                    kp = i
                    amax = absi
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k,k])
                for i = k+1:n
                    A[i,k] *= Akkinv
                end
            else
                throw(SingularException(k))
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:n
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    # ipiv2perm
    perm .= 1:n
    @inbounds begin
        for i in 1:n
            for j in 1:n
                B[i,j], B[ipiv[i],j] = B[ipiv[i],j], B[i,j]
            end
            perm[i], perm[ipiv[i]] = perm[ipiv[i]], perm[i]
        end
    end
    # Triangular solves
    @inbounds begin
        for k in 1:n
            for j in 1:n       # naivesub! for UnitLowerTriangular
                xj = B[j,k]
                for i in j+1:n
                    B[i,k] -= A[i,j] * xj
                end
            end
            for j in n:-1:1    # naivesub! for UpperTriangular
                xj = B[j,k] = A[j,j] \ B[j,k]
                for i in j-1:-1:1 # counterintuitively 1:j-1 performs slightly better
                    B[i,k] -= A[i,j] * xj
                end
            end
        end
    end
end

struct myexpmWorkspace{T}
    U::Matrix{T}
    V::Matrix{T}
    W::Matrix{T}
    A2::Matrix{T}
    A4::Matrix{T}
    A6::Matrix{T}
    ipiv::Vector{Int}
    perm::Vector{Int}
end
myexpmWorkspace(N, ::Type{T}) where {T} = myexpmWorkspace((zeros(T,N,N) for i in 1:6)..., zeros(Int,N), zeros(Int,N))

const myexpm_workspace_cache = Dict{Tuple{Int,DataType},Tuple{Matrix,myexpmWorkspace}}()
get_expm_workspace(N, ::Type{T}) where {T} = get!(myexpm_workspace_cache, (N,T)) do; (zeros(T,N,N), myexpmWorkspace(N,T)) end
function myexpm(A::StaticMatrix{N,N,T}) where {N,T}
    ws = get_expm_workspace(N,T)::Tuple{Matrix{T},myexpmWorkspace{T}}
    typeof(A)(myexpm!(copyto!(ws[1], A), ws[2]))
end
myexpm(A::Matrix{T}) where T<:BlasFloat = myexpm!(copy(A), myexpmWorkspace(size(A,1),T))
function myexpm!(A::Matrix{T}, ws::myexpmWorkspace{T}=myexpmWorkspace(size(A,1),T)) where T<:BlasFloat
    n = checksquare(A)
    nA = norm(A, 1)
    U, V, W, A2, A4, A6, ipiv, perm = ws.U, ws.V, ws.W, ws.A2, ws.A4, ws.A6, ws.ipiv, ws.perm
    ## For sufficiently small nA, use lower order Padé-Approximations
    if (nA <= 2.1)
        mul!(A2, A, A)
        if nA > 0.95
            @evalmatrixpoly!(V, W, 1, A2, T(8821612800), T(302702400), T(2162160), T(3960), T(1))
            mul!(U, A, V)
            @evalmatrixpoly!(V, W, 1, A2, T(17643225600), T(2075673600), T(30270240), T(110880), T(90))
        elseif nA > 0.25
            @evalmatrixpoly!(V, W, 1, A2, T(8648640), T(277200), T(1512), T(1))
            mul!(U, A, V)
            @evalmatrixpoly!(V, W, 1, A2, T(17297280), T(1995840), T(25200), T(56))
        elseif nA > 0.015
            @evalmatrixpoly!(V, W, 1, A2, T(15120), T(420), T(1))
            mul!(U, A, V)
            @evalmatrixpoly!(V, W, 1, A2, T(30240), T(3360), T(30))
        else
            @evalmatrixpoly!(V, W, 1, A2, T(60), T(1))
            mul!(U, A, V)
            @evalmatrixpoly!(V, W, 1, A2, T(120), T(12))
        end
        W .= V .- U
        A .= V .+ U
        mygesv!(W, A, ipiv, perm)
    else
        s  = log2(nA/5.4)   # power of 2 later reversed by squaring
        if s > 0
            si = ceil(Int,s)
            A ./= convert(T,2^si)
        end

        mul!(A2, A, A)
        mul!(A4, A2, A2)
        mul!(A6, A2, A4)

        @evalmatrixpoly!(V, W, 3, A2, A4, A6, T(32382376266240000), T(1187353796428800), T(10559470521600), T(33522128640), T(40840800), T(16380), T(1))
        mul!(U, A, V)
        @evalmatrixpoly!(V, W, 3, A2, A4, A6, T(64764752532480000), T(7771770303897600), T(129060195264000), T(670442572800), T(1323241920), T(960960), T(182))
        W .= V .- U
        A .= V .+ U
        mygesv!(W, A, ipiv, perm)

        if s > 0            # squaring to reverse dividing by power of 2
            for t=1:si
                mul!(W, A, A)
                copyto!(A, W)
            end
        end
    end
    A
end

# for T in (Complex{Float32}, Complex{Float64}, Float32, Float64)
#     A = rand(T,20,20)
#     nA = norm(A,1)
#     println((norm(expm(A) - myexpm(A)),
#              norm(expm(A.*1.5./nA) - myexpm(A.*1.5./nA)),
#              norm(expm(A.*0.75./nA) - myexpm(A.*0.75./nA)),
#              norm(expm(A.*0.1./nA) - myexpm(A.*0.1./nA)),
#              norm(expm(A.*0.01./nA) - myexpm(A.*0.01./nA))))
# end
