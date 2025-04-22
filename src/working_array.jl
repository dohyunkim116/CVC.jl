struct working_array{TT <: Real, T <: Real}
    xtz     :: Vector{T}
    xtz̃     :: Vector{T}
    xtz2    :: Vector{T}
    xtz̃2    :: Vector{T}
    xtystar :: Vector{T}
    xtystar2 :: Vector{T}
    v       :: Vector{T}
    v2      :: Vector{T}
    XkjtZ   :: Matrix{T}
    XkjtZ̃   :: Matrix{T}
    XkjtZ2  :: Matrix{T}
    XkjtZ̃2  :: Matrix{T}
    A       :: Matrix{T}
    A2      :: Matrix{T}
    Ã       :: Matrix{T}
    Ã2      :: Matrix{T}
    B       :: Matrix{T}
    B2      :: Matrix{T}
    B̃       :: Matrix{T}
    B̃2      :: Matrix{T}
    E       :: Matrix{T}
    E2      :: Matrix{T}
    Xₖʲ     :: Matrix{T}
    Xₖʲ2    :: Matrix{T}
    Xₖʲ_sla_array :: Array{SnpLinAlg{T}, 1}
    Z_kj    :: Matrix{T}
    Z̃_kj    :: Matrix{T}
    Q_kj    :: Vector{T}
end

function working_array(
    m                :: cvc{TT, T},
    k                :: Int
) where {TT, T}
n, q = size(m.H)
xtz       = Vector{T}(undef, m.M[k,1]) 
xtz̃       = similar(xtz) 
xtystar  = similar(xtz)
xtz2      = Vector{T}(undef, m.M[k,m.J]) 
xtz̃2      = similar(xtz2) 
xtystar2 = similar(xtz2)
v = similar(xtz)
v2 = similar(xtz2)
XkjtZ    = Matrix{T}(undef, m.M[k,1], m.B)
XkjtZ̃    = similar(XkjtZ)
XkjtZ2   = Matrix{T}(undef, m.M[k,m.J], m.B)
XkjtZ̃2   = similar(XkjtZ2)
A = Matrix{T}(undef, m.M[k,1], q)
A2 = Matrix{T}(undef, m.M[k,m.J], q)
Ã = similar(A)
Ã2 = similar(A2)
B = Matrix{T}(undef, m.M[k,1], n)
B2 = Matrix{T}(undef, m.M[k,m.J], n) 
B̃ = similar(B)
B̃2 = similar(B2)
E = Matrix{T}(undef, n, m.M[k,1])
E2 = Matrix{T}(undef, n, m.M[k,m.J])
Xₖʲ = Matrix{T}(undef, n, m.M[k,1])
Xₖʲ2 = Matrix{T}(undef, n, m.M[k,m.J])
Xₖʲ_sla_array = Array{SnpLinAlg{T}, 1}(undef, 1)
Z_kj = Matrix{T}(undef, n, m.B)
Z̃_kj = similar(Z_kj)
Q_kj = Vector{T}(undef, n)
working_array{TT, T}(
    xtz, xtz̃, xtz2, xtz̃2, xtystar, xtystar2, v, v2,
    XkjtZ, XkjtZ̃, XkjtZ2, XkjtZ̃2, A, A2, Ã, Ã2, 
    B, B2, B̃, B̃2, E, E2,
    Xₖʲ, Xₖʲ2, Xₖʲ_sla_array, Z_kj, Z̃_kj, Q_kj
)
end

function mmap_U_kj_array!(
    U_kj_array :: Array{AbstractArray{T}, 1},
    n::Int,
    B::Int,
    K::Int,
    temp_dir_path :: AbstractString;
    parallel :: Bool = false
) where T
    J = length(U_kj_array)
    if parallel
        Threads.@threads for j in 1:J
            io_U_jk = open("$(temp_dir_path)/U_jk_$(j).bin", "w+")
            U_kj_array[j] = Mmap.mmap(io_U_jk, Array{T,3}, (n, B, K))
            fill!(U_kj_array[j], zero(T))
        end
    else
        for j in 1:J
            io_U_jk = open("$(temp_dir_path)/U_jk_$(j).bin", "w+")
            U_kj_array[j] = Mmap.mmap(io_U_jk, Array{T,3}, (n, B, K))
            fill!(U_kj_array[j], zero(T))
        end
    end
end

function mmap_Ũ_kj_array!(
    Ũ_kj_array :: Array{AbstractArray{T}, 1},
    n::Int,
    B::Int,
    K::Int,
    temp_dir_path :: AbstractString;
    parallel :: Bool = false
) where T
    J = length(Ũ_kj_array)
    if parallel
        Threads.@threads for j in 1:J
            io_Ũ_jk = open("$(temp_dir_path)/Ũ_jk_$(j).bin", "w+")
            Ũ_kj_array[j] = Mmap.mmap(io_Ũ_jk, Array{T,3}, (n, B, K))
            fill!(Ũ_kj_array[j], zero(T))
        end
    else
        for j in 1:J
            io_Ũ_jk = open("$(temp_dir_path)/Ũ_jk_$(j).bin", "w+")
            Ũ_kj_array[j] = Mmap.mmap(io_Ũ_jk, Array{T,3}, (n, B, K))
            fill!(Ũ_kj_array[j], zero(T))
        end
    end
end

struct solve_working_array{TT <: Real, T <: Real}
    gramMat_array::Vector{Matrix{T}}
    rvec_array::Vector{Vector{T}}
    UmŨ_k_array::Vector{Matrix{T}}
    UmHHtU_l_array::Vector{Matrix{T}}
    HtU_jl_array::Vector{Matrix{T}}
    HHtU_jl_array::Vector{Matrix{T}}
    trKₖYstar_array::Vector{Vector{T}}
    trKₖPYstar_array::Vector{Vector{T}}
    trPKₖPYstar_array::Vector{Vector{T}}
end

function solve_working_array(
    m::cvc{TT, T}
) where {TT, T}
    gramMat_array = Vector{Matrix{T}}(undef, m.J + 1)
    rvec_array = Vector{Vector{T}}(undef, m.J + 1)
    UmŨ_k_array = similar(gramMat_array)
    UmHHtU_l_array = similar(gramMat_array)
    HtU_jl_array = similar(gramMat_array)
    HHtU_jl_array = similar(gramMat_array)
    trKₖYstar_array = similar(rvec_array)
    trKₖPYstar_array = similar(rvec_array)
    trPKₖPYstar_array = similar(rvec_array)
    n = size(m.H, 1)
    iszero(m.H) ? r = 0 : r = size(m.H, 2)
    for j in 1:(m.J + 1)
        gramMat_array[j] = Matrix{T}(undef, m.K + 1, m.K + 1)
        gramMat_array[j][m.K + 1, m.K + 1] = n - r
        UmŨ_k_array[j] = Matrix{T}(undef, n, m.B)
        UmHHtU_l_array[j] = Matrix{T}(undef, n, m.B)
        HtU_jl_array[j] = Matrix{T}(undef, size(m.H, 2), m.B)
        HHtU_jl_array[j] = Matrix{T}(undef, n, m.B)
        rvec_array[j]    = Vector{T}(undef, m.K + 1)
        trKₖYstar_array[j] = Vector{T}(undef, m.K)
        trKₖPYstar_array[j] = Vector{T}(undef, m.K)
        trPKₖPYstar_array[j] = Vector{T}(undef, m.K)
    end
    solve_working_array{TT, T}(
        gramMat_array, rvec_array, UmŨ_k_array, UmHHtU_l_array,
        HtU_jl_array, HHtU_jl_array, 
        trKₖYstar_array, trKₖPYstar_array, trPKₖPYstar_array
    )
end
