# This code has to be submitted for pull request to SnpArrays.jl 
# since the original code had an error with the imputation.
# function Base.getindex(s::SnpLinAlg{T}, i::Int, j::Int) where T
#    x = SnpArrays.convert(T, getindex(s.s, i, j), s.model)
#    s.impute && isnan(x) && (x = s.μ[j])
#    s.center && (x -= s.μ[j])
#    s.scale && (x *= s.σinv[j])
#    return x
# end

# mapsnp! functions are authored by Brendon Chau
"""
    mapsnp!(f, g, μ, σinv)

In-place additive imputation of a compressed column vector of a SnpArray `g` 
into vector of floating point numbers `f`. Centering and scaling are performed
automatically using scalars `μ` and `σinv`.
"""
@inline function mapsnp!(
    f::AbstractVector{T}, 
    g::AbstractVector{UInt8}, 
    μ::T, 
    σinv::T
    ) where {T<:AbstractFloat}
    K, rem = divrem(length(f), 4)
    @inbounds for j in 1:K
        # Base.Cartesian.@nexprs 4 i -> begin
        #     snp_i = (g[j] >> ((i - 1) * 2)) & 0x03
        #     f[(j - 1) << 2 + i] = (T((snp_i >= 0x02) * (snp_i - 0x01)) - !isone(snp_i) * μ) * σinv
        # end
        @simd for i in 1:4
            # rotate the correct value into the 0 and 1 place
            snp = (g[j] >> ((i - 1) << 1)) & 0x03
            # additive model imputation in-place
            f[(j - 1) << 2 + i] = (T((snp >= 0x02) * (snp - 0x01)) - !isone(snp) * μ) * σinv
        end
    end
    @inbounds if rem > 0
        @simd for i in 1:rem
            snp = (g[K + 1] >> ((i - 1) << 1)) & 0x03
            # additive model imputation in-place
            f[K << 2 + i] = (T((snp >= 0x02) * (snp - 0x01)) - !isone(snp) * μ) * σinv
        end
    end
end

"""
    mapsnp!(F, G, μ, σinv)

In-place additive imputation of a compressed SnpArray `G` into a matrix of 
floating point numbers `F`, parallelized over the columns of `F` and `G`. 
Centering and scaling are performed automatically using vectors `μ` and `σinv`.
"""
function mapsnp!(
    F::AbstractMatrix{T}, 
    G::AbstractMatrix{UInt8}, 
    μ::AbstractVector{T}, 
    σinv::AbstractVector{T}
    ) where {T<:AbstractFloat}
    @inbounds Threads.@threads for j in axes(F, 2)
        mapsnp!(view(F, :, j), view(G, :, j), μ[j], σinv[j])
    end
    return F
end

function mapsnp!(
    F::AbstractMatrix{T}, 
    G::AbstractSnpArray, 
    μ::AbstractVector{T}, 
    σinv::AbstractVector{T}
    ) where {T<:AbstractFloat}
    mapsnp!(F, G.data, μ, σinv)
    return F
end

function standardize!(X::SnpLinAlg, G::SnpArray)
    n, p = size(G)
    @inbounds @views for j in 1:p
        #X.μ[j] = (X.μ[j] * G.columncounts[2,j] + G.columncounts[3,j] + 2G.columncounts[4,j]) / n
        X.σinv[j] = (abs2(X.μ[j]) * G.columncounts[1, j] + 
            abs2(1.0 - X.μ[j]) * G.columncounts[3, j] + 
            abs2(2.0 - X.μ[j]) * G.columncounts[4,j]) / n
        X.σinv[j] = X.σinv[j] > 0 ? inv(sqrt(X.σinv[j])) : one(eltype(X.σinv))
    end
end

function standardizeW!(W)
    dt = StatsBase.fit(ZScoreTransform, W, dims = 1)
    StatsBase.transform!(dt, W)
    W[:,1] .= ones(eltype(W))
end