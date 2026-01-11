module CVC

using Distributions, LinearAlgebra, Printf, Reexport, Roots, SnpArrays, Random, StatsBase
using NonNegLeastSquares
using CodecZlib, CodecXz, CodecBzip2, CodecZstd, TranscodingStreams, Mmap
using DelimitedFiles: readdlm, writedlm
using DataFrames
import Base: length, Slice, OneTo, UnitRange
import LinearAlgebra: BlasReal, copytri!
import Statistics: mean, var
@reexport using StatsModels

include("km.jl")
include("snparray_tools.jl")
include("codec.jl")
include("synvar.jl")

export 
    # types
    cvc,
    working_array,
    KaplanMeier,
    # functions
    fit_me!,
    length,
    mean,
    nobs,
    Survival,
    var,
    standardize!,
    generate_design_matrices,
    generate_design_matrices_with_missing,
    generate_MAFs,
    compute_mean_component,
    parse_datapath,
    get_N,
    get_Mₖ,
    make_complete_data,
    coef,
    coefnames,
    show_sim_res,
    show_sim_res2,
    fit_res,
    save_sim_res,
    paramnames,
    params

function compute_M!(M::Matrix, M_k::Vector, K, J)
    for k in 1:K
        # All blocks except the last get the same size
        base_size = M_k[k] ÷ J
        
        # First J-1 blocks all get exactly base_size
        for j in 1:(J-1)
            M[k,j] = base_size
        end
        
        # Last block gets all remaining SNPs
        M[k,J] = M_k[k] - base_size * (J-1)
    end
end

function construct_sa_array!(sa_array::Array{SnpArray,1}, genopath::AbstractString, n::Int, K::Int)
    for k in 1:K
        Gpath = isfile("$(genopath)/G$k.bed.gz") ? "$(genopath)/G$k.bed.gz" : "$(genopath)/G$k.bed"
        sa_array[k] = SnpArray(Gpath, n)
    end
end

function filter_sa_array!(
    sa_array::Array{SnpArray,1},
    rowinds::AbstractVector{<:Integer}
    )
    n = size(sa_array[1], 1)
    if eltype(rowinds) == Bool
        @assert length(rowinds) == n "Length of rowinds should be same as number of individuals"
        rmask = rowinds
    else
        rmask = falses(n)
        rmask[rowinds] .= true
    end
    for k in eachindex(sa_array)
        Gₖ = SnpArray(undef, count(rmask), size(sa_array[k], 2))
        Gₖ .= @view sa_array[k][rmask, :]
        sa_array[k] = Gₖ
    end
end

function construct_sla_array!(sla_array::Array{SnpLinAlg{T}, 1}, sa_array::Array{SnpArray,1}, K::Int) where T
    for k in 1:K
        Xₖ = SnpLinAlg{T}(sa_array[k], model = ADDITIVE_MODEL, center = true, scale = true)
        standardize!(Xₖ, sa_array[k])
        sla_array[k] = Xₖ
    end
end

function construct_sla_array!(sla_array::Array{SnpLinAlg{T}, 1}, genopath, n::Int, K::Int) where T
    for k in 1:K
        Gpath = isfile("$(genopath)/G$k.bed.gz") ? "$(genopath)/G$k.bed.gz" : "$(genopath)/G$k.bed"
        Gₖ = SnpArray(Gpath, n)
        Xₖ = SnpLinAlg{T}(Gₖ, model = ADDITIVE_MODEL, center = true, scale = true)
        standardize!(Xₖ,Gₖ)
        sla_array[k] = Xₖ
    end
end

function construct_randZ̃!(
    randZ̃::Matrix{T}, 
    H::Union{Matrix{T}, 
            SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true},
            Transpose{T, Matrix{T}},
            Transpose{T, SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true}}
            },
    randZ::Matrix{T}
    ) where T    
    HtZ = Matrix{T}(undef, size(H, 2), size(randZ, 2))
    mul!(HtZ, transpose(H), randZ)
    mul!(randZ̃, H, HtZ)
end

function construct_randZ̃̃!(
    randZ̃̃::Matrix{T}, 
    H̃::Union{Matrix{T}, 
            SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true},
            Transpose{T, Matrix{T}},
            Transpose{T, SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true}}
            },
    H::Union{Matrix{T}, 
            Transpose{T, Matrix{T}},
            SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true},
            Transpose{T, SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true}}
            },
    randZ::Matrix{T}
    ) where T
    H̃tZ = Matrix{T}(undef, size(H̃, 2), size(randZ, 2))
    mul!(H̃tZ, transpose(H̃), randZ)
    mul!(randZ̃̃, H, H̃tZ)
end

# HHtystar1 is used in computing trYstarV, V_kj and trKₖPYstar
function update_HHtystar!(
    HHtystar1::Vector{T},
    H::Union{Matrix{T}, 
            SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true},
            Transpose{T, Matrix{T}},
            Transpose{T, SubArray{T, 2, Matrix{T}, Tuple{Slice{OneTo{Int64}}, UnitRange{Int64}}, true}}
            },
    ystar1::Vector{T}
    ) where T
    Htystar1 = Vector{T}(undef, size(H, 2))
    mul!(Htystar1, transpose(H), ystar1)    
    mul!(HHtystar1, H, Htystar1)
end

struct cvc{TT <: Real, T <: Real}
    # data
    ỹ           :: Vector{TT}    # observed event times: zᵢ = yᵢ ∧ cᵢ
    δ           :: Vector{Bool}  # censoring status: δᵢ = 1{yᵢ ≤ cᵢ}
    cr          :: Real
    model       :: Symbol
    datapath    :: String # data path containing plink1 genotype file set, and covaraite matrix file
    temp_dir    :: String # temporary directory for storing memory-mapped files
    N           :: Int
    Mtotal      :: Int
    C           :: Int
    K           :: Int
    J           :: Int
    B           :: Int
    M           :: Matrix{Int}
    M̃           :: Matrix{Int}
    M_k         :: Vector{Int}
    maf_mat_jack:: Matrix{T}
    # parameters
    isfitted    :: Vector{Bool}
    ϕg          :: Vector{TT}
    ϕe          :: Vector{TT}
    ϕ           :: Vector{TT}
    ϕse         :: Vector{TT}
    h2          :: Vector{TT}
    h2se        :: Vector{TT}
    h2k         :: Vector{TT}
    h2kse       :: Vector{TT}
    Jack        :: Matrix{T} # matrix containing jackknife estimates of ϕg and ϕe
    ϕg_nn       :: Vector{TT}
    ϕe_nn       :: Vector{TT}
    ϕ_nn        :: Vector{TT}
    ϕse_nn      :: Vector{TT}
    h2_nn       :: Vector{TT}
    h2se_nn     :: Vector{TT}
    h2k_nn      :: Vector{TT}
    h2kse_nn    :: Vector{TT}
    Jack_nn        :: Matrix{T}
    # synthetic variables
    ystar1      :: Vector{TT}      # first moment synthetic variable E(y⋆)=E(y)
    ystar2      :: Vector{TT}      # second moment synthetic variable E(y⋆²)=E(y²)
    # KM objects
    km          :: KaplanMeier{TT, T} # Kaplan-Meier estiamte of phenotype distn
    kmc         :: KaplanMeier{TT, T} # Kaplan-Meier estiamte of censoring distn
    cec         :: Vector{TT}      # compensation term ∫ᵗʲ G(t) / (1 - G(t)) dt in first moment Leurgan's method
    cec2        :: Vector{TT}      # compensation term ∫ᵗʲ 2tG(t) / (1 - G(t)) dt in second moment Leurgan's method
    # fixed effect related variables (these are not used in the current version of the code)
    ŷ           :: Vector{TT}
    w           :: Matrix{TT}
    # for solving normal equations
    trYstarV    :: Vector{T}
    # working arrays
    H           :: Matrix{T}
    HHtystar1   :: Vector{T}
    H̃           :: Matrix{T}
    d           :: Vector{T}
    randZ       :: Matrix{T}
    randZ̃       :: Matrix{T}
    #randZ̃̃       :: Matrix{T}
    L̃_0        :: Vector{T}
    L̃        :: Matrix{T}
    R̃_0        :: Vector{T}
    R̃        :: Matrix{T}
    S̃_0        :: Vector{T}
    S̃        :: Matrix{T}
    F̃_0        :: Vector{T}
    F̃        :: Matrix{T}
    Ṽ_0        :: Vector{T}
    Ṽ        :: Matrix{T}
    Q̃_0        :: Matrix{T}
    Q̃_kj_array  :: Array{AbstractArray{T},1}
    U_0        :: Array{T,3}
    Ũ_0        :: Array{T,3}
    HHtU_0     :: Array{T,3}
    U_kj_array :: Array{AbstractArray{T},1}
    Ũ_kj_array :: Array{AbstractArray{T},1}
    HtU_jl       :: Matrix{T}
    HHtU_jl     :: Matrix{T}
    ŨZ̃̃          :: Array{T,3}
    trPKₖPYstar_array :: Array{T, 3}
    sa_array   :: Array{SnpArray, 1}
    sla_array   :: Array{SnpLinAlg, 1}
end

function cvc(
    ỹ            :: AbstractVector{TT},
    δ            :: Union{Vector{Bool}, BitVector},
    w            :: AbstractArray,
    datapath     :: String,
    tempdir_parent     :: AbstractString;
    model = :tte,
    J = 100,     # Allow customizing initial number of jackknife blocks
    B = 10       # Allow customizing number of random vectors
    ) where {TT}
    @assert length(ỹ) == length(δ) "Lengths of ỹ and δ should be same"
    T        = Float32
    N        = get_N(datapath)
    K        = filter(contains(r"G\d+.bed"), readdir(datapath)) |> length
    sa_array = Array{SnpArray,1}(undef, K)
    construct_sa_array!(sa_array, datapath, N, K)
    censored = true
    cr = round((N - count(δ)) / N, digits = 4)
    sla_array = Array{SnpLinAlg{T}, 1}(undef, K);
    construct_sla_array!(sla_array, sa_array, K)
    M_k      = Vector{Int}(undef, K)
    M_k     .= get_Mₖ.(1:K, datapath)
    Mtotal = sum(M_k)
    C = size(w, 2)
    # Dynamically adjust J if any component has fewer SNPs than J
    min_M_k = minimum(M_k)
    @assert min_M_k >= 2 "Each component must have at least 2 SNPs"
    if min_M_k <= J
        J_old = J
        J = max(2, min_M_k)  # Ensure J doesn't go below 2
        @info "Reducing jackknife blocks from $J_old to $J because smallest component has only $min_M_k SNPs"
    end
    
    M        = Matrix{Int}(undef, K, J)
    compute_M!(M, M_k, K, J)
    M̃        = similar(M)
    M̃       .= M_k .- M
    maf_mat_jack = Matrix{T}(undef, K, J + 1)
    ŷ        = Vector{TT}(undef, N)
    if !iszero(w)
        w = convert(Matrix{Float64}, w)
        F        = svd(w)
        r         = rank(Diagonal(F.S))
        H         = @view(F.U[:, 1:r])
    else
        H = w
    end
    H         = convert(Matrix{T}, H)
    temp_dir_path = mktempdir(tempdir_parent, cleanup = false)
    isfitted  = [false]
    ϕg        = Vector{T}(undef, K)
    ϕe        = Vector{T}(undef, 1)
    ϕ         = Vector{T}(undef, K + 1)
    ϕse       = Vector{T}(undef, K + 1)
    h2        = Vector{T}(undef, 1)
    h2se      = Vector{T}(undef, 1)
    h2k       = Vector{T}(undef, K + 1)
    h2kse     = Vector{T}(undef, K + 1)
    Jack      = Matrix{T}(undef, K + 1, J + 1)
    ϕg_nn     = Vector{T}(undef, K)
    ϕe_nn     = Vector{T}(undef, 1)
    ϕ_nn      = Vector{T}(undef, K + 1)
    ϕse_nn    = Vector{T}(undef, K + 1)
    h2_nn     = Vector{T}(undef, 1)
    h2se_nn   = Vector{T}(undef, 1)
    h2k_nn       = Vector{T}(undef, K + 1)
    h2kse_nn     = Vector{T}(undef, K + 1)
    Jack_nn   = Matrix{T}(undef, K + 1, J + 1)
    ystar1    = Vector{T}(undef, N)
    ystar2    = Vector{T}(undef, N)
    km        = KaplanMeier(ỹ, δ, TS = T)
    kmc       = KaplanMeier(ỹ, .!δ, TS = T)
    cec = intgrl(kmc)
    cec2 = intgrl2(kmc)
    trYstarV  = Vector{T}(undef, 1)
    HHtystar1 = Vector{T}(undef, N)
    H̃         = similar(H)
    d         = Vector{T}(undef, N)
    randZ    = Matrix{T}(undef, N, B)
    randZ̃    = Matrix{T}(undef, N, B)
    L̃_0     = zeros(T, K)
    L̃     = zeros(T, K, J)
    R̃_0     = zeros(T, K)
    R̃     = zeros(T, K, J)
    S̃_0     = zeros(T, K)
    S̃     = zeros(T, K, J)
    F̃_0     = zeros(T, K)
    F̃     = zeros(T, K, J)
    Ṽ_0     = zeros(T, K)
    Ṽ     = zeros(T, K, J)
    Q̃_0     = Matrix{T}(undef, N, K)
    fill!(Q̃_0, zero(T))
    Q̃_kj_array = Array{Matrix{T},1}(undef, J)
    for j in 1:J
        Q̃_kj_array[j] = Matrix{T}(undef, N, K)
        fill!(Q̃_kj_array[j], zero(T))
    end
    
    U_0 = Array{T,3}(undef, N, B, K)
    fill!(U_0, zero(T))
    Ũ_0 = Array{T,3}(undef, N, B, K)
    fill!(Ũ_0, zero(T))
    HHtU_0 = Array{T,3}(undef, N, B, K)
    fill!(HHtU_0, zero(T))

    U_kj_array = Array{AbstractArray{T}, 1}(undef, J)
    Ũ_kj_array = Array{AbstractArray{T}, 1}(undef, J)
    mmap_U_kj_array!(U_kj_array, N, B, K, temp_dir_path)
    mmap_Ũ_kj_array!(Ũ_kj_array, N, B, K, temp_dir_path)

    HtU_jl    = Matrix{T}(undef, size(H, 2), B)
    HHtU_jl  = similar(randZ)

    ŨZ̃̃ = similar(Ũ_0)
    trPKₖPYstar_array = Array{T, 3}(undef, 1, 1, K)
    synthetic!(ystar1, ystar2, ỹ, δ, kmc, cec, cec2)
    ystar1 = convert(Vector{T}, ystar1)
    ystar2 = convert(Vector{T}, ystar2)
    d .= ystar2 .- ystar1.^2
    H̃ .= d .* H
    rand!(Normal(), randZ)
    construct_randZ̃!(randZ̃, H, randZ)
    update_HHtystar!(HHtystar1, H, ystar1)

    cvc{TT, T}(
        ỹ, δ, cr, model, datapath, temp_dir_path, N, Mtotal, C, K, J, B, M, M̃, M_k, maf_mat_jack,
        isfitted, ϕg, ϕe, ϕ, ϕse, h2, h2se, h2k, h2kse, Jack,
        ϕg_nn, ϕe_nn, ϕ_nn, ϕse_nn, h2_nn, h2se_nn, h2k_nn, h2kse_nn, Jack_nn,
        ystar1, ystar2, km, kmc, cec, cec2,
        ŷ, w, trYstarV,
        H, HHtystar1, H̃, d, randZ, randZ̃,
        L̃_0, L̃, R̃_0, R̃, S̃_0, S̃, F̃_0, F̃, Ṽ_0, Ṽ, 
        Q̃_0, Q̃_kj_array, 
        U_0, Ũ_0, HHtU_0,
        U_kj_array, Ũ_kj_array,
        HtU_jl, HHtU_jl, ŨZ̃̃, trPKₖPYstar_array,
        sa_array, sla_array
        )
end

function params_pt(m::cvc; nnls = false)
    if nnls
        pt_est_vec = [m.h2k_nn; m.h2_nn]
    else 
        pt_est_vec = [m.h2k; m.h2]
    end
    pt_est_vec
end

function params_se(m::cvc; nnls = false)
    if nnls
        se_est_vec = [m.h2kse_nn; m.h2se_nn]
    else 
        se_est_vec = [m.h2kse; m.h2se]
    end
    se_est_vec
end

function params(m::cvc; nnls = false)
    [params_pt(m; nnls = nnls); params_se(m; nnls = nnls)]
end

paramnames_pt(m::cvc) = [["h2_$k" for k in 1:m.K]; "h2_e"; "h2"]
paramnames_se(m::cvc) = [["h2_$(k)_se" for k in 1:m.K]; "h2_e_se"; "h2_se"]
paramnames(m::cvc) = [paramnames_pt(m); paramnames_se(m)]

function paramtable(m::cvc; nnls = false)
    if nnls
        pars = params(m; nnls)
    else
        pars = params(m; nnls)
    end
    hcat(paramnames(m), pars)
end

coefnames(m::cvc)= [["h2_$k" for k in 1:m.K]; "h2_e"; "h2"]

function coef(m::cvc; nnls = false)
    if nnls
        hcat([m.h2k_nn; m.h2_nn],[m.h2kse_nn; m.h2se_nn])
    else
        hcat([m.h2k; m.h2],[m.h2kse; m.h2se])
    end
end
nobs(m::cvc) = length(m.ỹ)
ncens(m::cvc) = nobs(m) - count(m.δ)
propcens(m::cvc) = m.model == :tte ? round(ncens(m) / nobs(m), sigdigits = 2) : 0
ngpred(m::cvc) = sum(m.M_k)
nfpred(m::cvc) = iszero(m.w) ? 0 : size(m.w, 2)
ngcomp(m::cvc) = m.K

function coeftable(m::cvc; nnls = false)
    if nnls
        mcoefs = coef(m; nnls)
    else
        mcoefs = coef(m)
    end
    StatsModels.CoefTable(
        mcoefs,
        ["Estimate", "SE"],
        coefnames(m)
        )
end

function Base.show(io::IO, m::cvc)
    if !m.isfitted[1]
        @warn("The model has not been fit.")
        return nothing
    end
    println(io)
    println(io, "Censored variance component model by synthetic variables")
    println(io)
    println(io, "Number of observations: $(nobs(m))")
    println(io, "Right-censored: $(ncens(m)) ($(round(propcens(m) * 100, digits = 2))% of observations)")
    println(io, "Number of genotypes: $(ngpred(m))")
    println(io, "Number of covariates: $(nfpred(m)) (including the intercept)")
    println(io, "Number of genetic components: $(ngcomp(m))")
    println(io)
    println(io, "Regression coefficients without nonnegativity constraints:")
    show(io, coeftable(m))
    println(io)
    println(io)
    println(io, "Regression coefficients with nonnegativity constraints:")
    show(io, coeftable(m; nnls = true))
    println(io)
end

function get_start_end_Gₖʲ(m::cvc, k, j)
    s = j == m.J ? m.M_k[k] - m.M[k,j] + 1 : 1 + (j - 1) * m.M[k,j]
    e = s + m.M[k,j] - 1
    s, e
end

function update_maf_mat_jack!(m::cvc{TT, T}) where {TT, T}
    mafk_array = Array{Vector{T}, 1}(undef, m.K)
    for k in 1:m.K
        mafk_array[k] = maf(m.sa_array[k])
    end
    maf_sum_k = sum.(mafk_array)
    temp_mat = Matrix{T}(undef, m.J, m.K)
    for k in 1:m.K
        for j in 1:m.J
            s, e = get_start_end_Gₖʲ(m, k, j)
            temp_mat[j,k] = sum(@view(mafk_array[k][s:e]))
        end
    end
    temp_mat .= (Transpose(maf_sum_k) .- temp_mat) ./ Transpose(m.M̃)
    copyto!(m.maf_mat_jack, Transpose(temp_mat))
    @views m.maf_mat_jack[:, m.J + 1] .= maf_sum_k ./ m.M_k
    return m.maf_mat_jack
end

function __set_estimates!(
    m::cvc{TT, T},
    nnls::Bool,
    Jack::Matrix{T},
    h2_rowvec::Matrix{T},
    h2e_rowvec::Matrix{T},
    H2::Matrix{T},
    E::Matrix{T},
    B::Matrix{T}
    ) where {TT, T}
    if nnls
        m.ϕg_nn .= Jack[1:m.K, m.J + 1]
        m.ϕe_nn .= Jack[m.K+1, m.J + 1]
        m.ϕ_nn .= [m.ϕg_nn ; m.ϕe_nn]
        copyto!(m.ϕse_nn, std(@view(Jack[:, 1:m.J]), dims = 2, corrected = false) .* sqrt(m.J - 1))
        m.h2_nn .= h2_rowvec[1, m.J + 1]
        m.h2se_nn .= std(@view(h2_rowvec[1, 1:m.J]), corrected = false) .* sqrt(m.J - 1)
        copyto!(m.h2k_nn, @view(H2[:, m.J+1]))
        m.h2k_nn[m.K + 1] = h2e_rowvec[1, m.J + 1]
        copyto!(m.h2kse_nn, std(@view(H2[:,1:m.J]), dims = 2, corrected = false) .* sqrt(m.J - 1))
        m.h2kse_nn[m.K + 1] = std(@view(h2e_rowvec[1, 1:m.J]), corrected = false) * sqrt(m.J - 1)
    else
        m.ϕg .= Jack[1:m.K, m.J + 1]
        m.ϕe .= Jack[m.K+1, m.J + 1]
        m.ϕ .= [m.ϕg ; m.ϕe]
        copyto!(m.ϕse, std(@view(Jack[:, 1:m.J]), dims = 2, corrected = false) .* sqrt(m.J - 1))
        m.h2 .= h2_rowvec[1, m.J + 1]
        m.h2se .= std(@view(h2_rowvec[1, 1:m.J]), corrected = false) .* sqrt(m.J - 1)
        copyto!(m.h2k, @view(H2[:, m.J+1]))
        m.h2k[m.K + 1] = h2e_rowvec[1, m.J + 1]
        copyto!(m.h2kse, std(@view(H2[:,1:m.J]), dims = 2, corrected = false) .* sqrt(m.J - 1))
        m.h2kse[m.K + 1] = std(@view(h2e_rowvec[1, 1:m.J]), corrected = false) * sqrt(m.J - 1)
    end
end

function _set_estimates!(
    m::cvc{TT, T},
    nnls::Bool;
    overlap_annot::Bool = false,
    Ĩ::Matrix{T} = Matrix{T}(undef, sum(m.M_k), m.K)
    ) where {TT, T}
    if nnls
        Jack = m.Jack_nn
    else
        Jack = m.Jack
    end
    𝚽g = @view(Jack[1:m.K, :])
    totalgvar_rowvec = sum(𝚽g, dims = 1)
    totalvar_rowvec = sum(Jack, dims = 1)
    totalevar_rowvec = totalvar_rowvec .- totalgvar_rowvec
    h2_rowvec = Matrix{T}(undef, 1, m.J + 1)
    h2_rowvec .= totalgvar_rowvec ./ totalvar_rowvec
    h2e_rowvec = Matrix{T}(undef, 1, m.J + 1)
    h2e_rowvec .= totalevar_rowvec ./ totalvar_rowvec
    M̃plus = hcat(m.M̃, m.M_k)
    M_rowvec = sum(M̃plus, dims = 1)
    H2 = similar(𝚽g)
    if !overlap_annot
        H2 .= 𝚽g ./ totalvar_rowvec
    else
        p = sum(m.M_k)
        Ĩ𝚽g = Matrix{T}(undef, p, m.J + 1) # we want to import Ĩ as memory-mapped array
        mul!(Ĩ𝚽g, Ĩ, 𝚽g)
        𝚽tildeg = similar(𝚽g)
        𝚽tildeg .= 𝚽tildeg ./ M̃plus
        mul!(𝚽tildeg, transpose(Ĩ), Ĩ𝚽g)
        H2 .= 𝚽tildeg ./ totalvar_rowvec
    end
    E = H2 ./ h2_rowvec ./ M̃plus .* M_rowvec
    update_maf_mat_jack!(m)
    B = (H2 .* totalvar_rowvec) ./ M̃plus ./ (2 .* m.maf_mat_jack .* (1 .- m.maf_mat_jack))
    
    __set_estimates!(m, nnls, Jack, h2_rowvec, h2e_rowvec, H2, E, B)
end

function set_estimates!(m::cvc{TT, T}) where {TT, T}
    _set_estimates!(m, false)
    _set_estimates!(m, true)
end

include("fit.jl")
include("working_array.jl")
include("routine.jl")

end
