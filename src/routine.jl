function update_Xₖʲ!(
    m::cvc{TT, T}, 
    wa::working_array{TT, T}, 
    k::Int,
    j::Int,
    sla::Bool = false
    ) where {TT, T}
    s, e = get_start_end_Gₖʲ(m, k, j)
    if sla
        n = size(m.H, 1)
        Gₖʲ = SnpArray(undef, n, e - s + 1)
        Gₖʲ.data .= @view m.sa_array[k].data[:, s:e]
        wa.Xₖʲ_sla_array[1] = SnpLinAlg{T}(Gₖʲ, model = ADDITIVE_MODEL, center = true, scale = true)
        standardize!(wa.Xₖʲ_sla_array[1],Gₖʲ)
    else
        if j == m.J
            #copyto!(wa.Xₖʲ2, @view m.sla_array[k][:, s:e]);
            mapsnp!(wa.Xₖʲ2, @view(m.sa_array[k].data[:, s:e]), @view(m.sla_array[k].μ[s:e]), @view(m.sla_array[k].σinv[s:e]))
        else
            #copyto!(wa.Xₖʲ, @view m.sla_array[k][:, s:e]);
            mapsnp!(wa.Xₖʲ, @view(m.sa_array[k].data[:, s:e]), @view(m.sla_array[k].μ[s:e]), @view(m.sla_array[k].σinv[s:e]))
        end
    end
end
 
#wa.A is used in computing R_kj
function compute_L_kj(m::cvc{TT, T}, wa::working_array{TT, T}, j, sla::Bool = false) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            mul!(wa.A, transpose(wa.Xₖʲ_sla_array[1]), m.H)
            dot(wa.A, wa.A)
        elseif j == m.J
            mul!(wa.A2, transpose(wa.Xₖʲ_sla_array[1]), m.H)
            dot(wa.A2, wa.A2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            mul!(wa.A, transpose(wa.Xₖʲ), m.H)
            dot(wa.A, wa.A)
        elseif j == m.J
            mul!(wa.A2, transpose(wa.Xₖʲ2), m.H)
            dot(wa.A2, wa.A2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function compute_R_kj(m::cvc{TT, T}, wa::working_array{TT, T}, j, sla::Bool = false) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            mul!(wa.Ã, transpose(wa.Xₖʲ_sla_array[1]), m.H̃)
            dot(wa.A, wa.Ã)
        elseif j == m.J
            mul!(wa.Ã2, transpose(wa.Xₖʲ_sla_array[1]), m.H̃)
            dot(wa.A2, wa.Ã2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            mul!(wa.Ã, transpose(wa.Xₖʲ), m.H̃)
            dot(wa.A, wa.Ã)
        elseif j == m.J
            mul!(wa.Ã2, transpose(wa.Xₖʲ2), m.H̃)
            dot(wa.A2, wa.Ã2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function compute_S_kj(m::cvc{TT, T}, wa::working_array{TT, T}, j) where {TT, T}
    if j >= 1 && j < m.J
        mul!(wa.B, wa.A, transpose(m.H))
        mul!(wa.B̃, wa.A, transpose(m.H̃))
        dot(wa.B, wa.B̃)
    elseif j == m.J
        mul!(wa.B2, wa.A2, transpose(m.H))
        mul!(wa.B̃2, wa.A2, transpose(m.H̃))
        dot(wa.B2, wa.B̃2)
    else
        throw(ArgumentError("j must be between 1 and $(m.J)"))
    end
end

function compute_F_kj(m::cvc{TT, T}, wa::working_array{TT, T}, j, sla::Bool = false) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            wa.E .= m.d .* (wa.Xₖʲ_sla_array[1] .* wa.Xₖʲ_sla_array[1])
            sum(wa.E)
        elseif j == m.J
            wa.E2 .= m.d .* (wa.Xₖʲ_sla_array[1] .* wa.Xₖʲ_sla_array[1])
            sum(wa.E2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            wa.E .= m.d .* (wa.Xₖʲ .* wa.Xₖʲ)
            sum(wa.E)
        elseif j == m.J
            wa.E2 .= m.d .* (wa.Xₖʲ2 .* wa.Xₖʲ2)
            sum(wa.E2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function compute_V_kj(m::cvc{TT, T}, wa::working_array{TT, T}, j, sla::Bool = false) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            #m.HHtystar1 is used in update_trYstarV!()
            mul!(wa.v,transpose(wa.Xₖʲ_sla_array[1]),m.HHtystar1) 
            dot(wa.v,wa.v)
        elseif j == m.J
            mul!(wa.v2,transpose(wa.Xₖʲ_sla_array[1]),m.HHtystar1)
            dot(wa.v2,wa.v2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            #m.HHtystar1 is used in update_trYstarV!()
            mul!(wa.v,transpose(wa.Xₖʲ),m.HHtystar1) 
            dot(wa.v,wa.v)
        elseif j == m.J
            mul!(wa.v2,transpose(wa.Xₖʲ2),m.HHtystar1)
            dot(wa.v2,wa.v2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function update_Z_kj!(
    m::cvc{TT, T},
    wa::working_array{TT, T},
    j::Int,
    sla::Bool = false
    ) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            mul!(wa.XkjtZ, transpose(wa.Xₖʲ_sla_array[1]), m.randZ)
            mul!(wa.Z_kj, wa.Xₖʲ_sla_array[1], wa.XkjtZ)
        elseif j == m.J
            mul!(wa.XkjtZ2, transpose(wa.Xₖʲ_sla_array[1]), m.randZ)
            mul!(wa.Z_kj, wa.Xₖʲ_sla_array[1], wa.XkjtZ2)
        else 
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            mul!(wa.XkjtZ, transpose(wa.Xₖʲ), m.randZ)
            mul!(wa.Z_kj, wa.Xₖʲ, wa.XkjtZ)
        elseif j == m.J
            mul!(wa.XkjtZ2, transpose(wa.Xₖʲ2), m.randZ)
            mul!(wa.Z_kj, wa.Xₖʲ2, wa.XkjtZ2)
        else 
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function update_Z̃_kj!(
    m::cvc{TT, T},
    wa::working_array{TT, T},
    j::Int,
    sla::Bool = false
    ) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            mul!(wa.XkjtZ̃, transpose(wa.Xₖʲ_sla_array[1]), m.randZ̃)
            mul!(wa.Z̃_kj, wa.Xₖʲ_sla_array[1], wa.XkjtZ̃)
        elseif j == m.J
            mul!(wa.XkjtZ̃2, transpose(wa.Xₖʲ_sla_array[1]), m.randZ̃)
            mul!(wa.Z̃_kj, wa.Xₖʲ_sla_array[1], wa.XkjtZ̃2)
        else 
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            mul!(wa.XkjtZ̃, transpose(wa.Xₖʲ), m.randZ̃)
            mul!(wa.Z̃_kj, wa.Xₖʲ, wa.XkjtZ̃)
        elseif j == m.J
            mul!(wa.XkjtZ̃2, transpose(wa.Xₖʲ2), m.randZ̃)
            mul!(wa.Z̃_kj, wa.Xₖʲ2, wa.XkjtZ̃2)
        else 
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end
end

function update_Q_kj!(m::cvc{TT, T}, wa::working_array{TT, T}, j, sla::Bool = false) where {TT, T}
    if sla
        if j >= 1 && j < m.J
            mul!(wa.xtystar,transpose(wa.Xₖʲ_sla_array[1]),m.ystar1)
            mul!(wa.Q_kj,wa.Xₖʲ_sla_array[1],wa.xtystar)
        elseif j == m.J
            mul!(wa.xtystar2,transpose(wa.Xₖʲ_sla_array[1]),m.ystar1)
            mul!(wa.Q_kj,wa.Xₖʲ_sla_array[1],wa.xtystar2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    else
        if j >= 1 && j < m.J
            mul!(wa.xtystar,transpose(wa.Xₖʲ),m.ystar1)
            mul!(wa.Q_kj,wa.Xₖʲ,wa.xtystar)
        elseif j == m.J
            mul!(wa.xtystar2,transpose(wa.Xₖʲ2),m.ystar1)
            mul!(wa.Q_kj,wa.Xₖʲ2,wa.xtystar2)
        else
            throw(ArgumentError("j must be between 1 and $(m.J)"))
        end
    end 
end

function accumulate_L̃!(m::cvc{TT, T}, L_kj, k, j) where {TT, T}
    m.L̃_0[k] += L_kj
    m.L̃[k, j] -= L_kj
end

function accumulate_R̃!(m::cvc{TT, T}, R_kj, k, j) where {TT, T}
    m.R̃_0[k] += R_kj
    m.R̃[k, j] -= R_kj
end

function accumulate_S̃!(m::cvc{TT, T}, S_kj, k, j) where {TT, T}
    m.S̃_0[k] += S_kj
    m.S̃[k, j] -= S_kj
end

function accumulate_F̃!(m::cvc{TT, T}, F_kj, k, j) where {TT, T}
    m.F̃_0[k] += F_kj
    m.F̃[k, j] -= F_kj
end

function accumulate_Ṽ!(m::cvc{TT, T}, V_kj, k, j) where {TT, T}
    m.Ṽ_0[k] += V_kj
    m.Ṽ[k, j] -= V_kj
end

function accumulate_U!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j) where {TT, T}
    @views m.U_0[:,:,k] .+= wa.Z_kj
    @views m.U_kj_array[j][:,:,k] .-= wa.Z_kj
end

function accumulate_Ũ!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j) where {TT, T}
    @views m.Ũ_0[:,:,k] .+= wa.Z̃_kj
    @views m.Ũ_kj_array[j][:,:,k] .-= wa.Z̃_kj
end

function accumulate_Q̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j) where {TT, T}
    @views m.Q̃_0[:,k] += wa.Q_kj
    @views m.Q̃_kj_array[j][:,k] .-= wa.Q_kj
end

function update_L̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    L_kj = compute_L_kj(m, wa, j, sla)
    accumulate_L̃!(m, L_kj, k, j)
end

function update_R̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    R_kj = compute_R_kj(m, wa, j, sla)
    accumulate_R̃!(m, R_kj, k, j)
end

function update_S̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j) where {TT, T}
    S_kj = compute_S_kj(m, wa, j)
    accumulate_S̃!(m, S_kj, k, j)
end

function update_F̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    F_kj = compute_F_kj(m, wa, j, sla)
    accumulate_F̃!(m, F_kj, k, j)
end

function update_Ṽ!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    V_kj = compute_V_kj(m, wa, j, sla)
    accumulate_Ṽ!(m, V_kj, k, j)
end

function update_U!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    update_Z_kj!(m, wa, j, sla)
    accumulate_U!(m, wa, k, j)
end

function update_Ũ!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    update_Z̃_kj!(m, wa, j, sla)
    accumulate_Ũ!(m, wa, k, j)
end

function update_Q̃!(m::cvc{TT, T}, wa::working_array{TT, T}, k, j, sla) where {TT, T}
    update_Q_kj!(m, wa, j, sla)
    accumulate_Q̃!(m, wa, k, j)
end

function update_HHtU_k0!(m::cvc{TT, T}) where {TT, T}
    HtU_k0 = Matrix{T}(undef, size(m.H, 2), m.B)
    @views for k in 1:m.K
        mul!(HtU_k0, transpose(m.H), m.U_0[:,:,k])
        mul!(m.HHtU_0[:,:,k], m.H, HtU_k0)
    end
end

function update_T!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    if j == m.J + 1
        for l in 1:m.K
            @views swa.UmHHtU_l_array[j] .= m.U_0[:,:,l] .- m.HHtU_0[:,:,l]
            for k in l:m.K
                @views swa.UmŨ_k_array[j] .= m.U_0[:,:,k] .- m.Ũ_0[:,:,k]
                swa.gramMat_array[j][k,l] = dot(swa.UmŨ_k_array[j], swa.UmHHtU_l_array[j]) / m.B / m.M_k[k] / m.M_k[l]  
            end
        end
    elseif j >= 1 && j <= m.J
        for l in 1:m.K
            mul!(swa.HtU_jl_array[j], transpose(m.H), @view(m.U_kj_array[j][:,:,l]))
            mul!(swa.HHtU_jl_array[j], m.H, swa.HtU_jl_array[j])
            @views swa.UmHHtU_l_array[j] .= m.U_kj_array[j][:,:,l] .- swa.HHtU_jl_array[j]
            for k in l:m.K
                @views swa.UmŨ_k_array[j] .= m.U_kj_array[j][:,:,k] .- m.Ũ_kj_array[j][:,:,k]
                swa.gramMat_array[j][k,l] = dot(swa.UmŨ_k_array[j], swa.UmHHtU_l_array[j]) / m.B / m.M̃[k,j] / m.M̃[l,j]
            end
        end
    else
        throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
    end
end

function update_b!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    n = size(m.H, 1)
    if j == m.J + 1
        @views swa.gramMat_array[j][1:m.K, m.K + 1] .= n .- m.L̃_0 ./ m.M_k
    elseif j >= 1 && j <= m.J
        @views swa.gramMat_array[j][1:m.K, m.K + 1] .= n .- m.L̃[:,j] ./ m.M̃[:,j]
    else
        throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
    end
end

function update_trKₖYstar!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    if j == m.J + 1
        mul!(swa.trKₖYstar_array[j], transpose(m.Q̃_0), m.ystar1)
        swa.trKₖYstar_array[j] .= (swa.trKₖYstar_array[j] .+ m.F̃_0) ./ m.M_k
    elseif j >= 1 && j <= m.J
        mul!(swa.trKₖYstar_array[j], transpose(m.Q̃_kj_array[j]), m.ystar1)
        @views swa.trKₖYstar_array[j] .= (swa.trKₖYstar_array[j] .+ m.F̃[:,j]) ./ m.M̃[:,j]
    else
        throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
    end
end

function update_trKₖPYstar!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    if j == m.J + 1
        mul!(swa.trKₖPYstar_array[j], transpose(m.Q̃_0), m.HHtystar1)
        swa.trKₖPYstar_array[j] .= (swa.trKₖPYstar_array[j] .+ m.R̃_0) ./ m.M_k
    elseif j >= 1 && j <= m.J
        mul!(swa.trKₖPYstar_array[j], transpose(m.Q̃_kj_array[j]), m.HHtystar1)
        @views swa.trKₖPYstar_array[j] .= (swa.trKₖPYstar_array[j] .+ m.R̃[:,j]) ./ m.M̃[:,j]
    else
        throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
    end
end

# function update_trPKₖPYstar!(
#     m::cvc{TT, T},
#     j::Int
#     ) where {TT, T}
#     if j == m.J + 1
#         m.ŨZ̃̃ .= m.Ũ_0 .* m.randZ̃̃ ./ m.B
#         sum!(m.trPKₖPYstar_array, m.ŨZ̃̃)
#         copyto!(m.trPKₖPYstar, m.trPKₖPYstar_array)
#         m.trPKₖPYstar .= (m.trPKₖPYstar .+ m.Ṽ_0) ./ m.M_k
#     elseif j >= 1 && j <= m.J
#         @views m.ŨZ̃̃ .= m.Ũ_kj_array[j] .* m.randZ̃̃ ./ m.B
#         sum!(m.trPKₖPYstar_array, m.ŨZ̃̃)
#         copyto!(m.trPKₖPYstar, m.trPKₖPYstar_array)
#         @views m.trPKₖPYstar .= (m.trPKₖPYstar .+ m.Ṽ[:,j]) ./ m.M̃[:,j]
#     else
#         throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
#     end
# end

function new_update_trPKₖPYstar!(
    m::cvc{TT, T},
    swa::solve_working_array,
    j::Int
    ) where {TT, T}
    if j == m.J + 1
        swa.trPKₖPYstar_array[j] .= (m.S̃_0 .+ m.Ṽ_0) ./ m.M_k
    elseif j >= 1 && j <= m.J
        @views swa.trPKₖPYstar_array[j] .= (m.S̃[:,j] .+ m.Ṽ[:,j]) ./ m.M̃[:,j]
    else
        throw(ArgumentError("j must be between 1 and $(m.J + 1)"))
    end
end

function update_c!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    update_trKₖYstar!(m, swa, j)
    update_trKₖPYstar!(m, swa, j)
    #update_trPKₖPYstar!(m, j)
    new_update_trPKₖPYstar!(m, swa, j)
    @views swa.rvec_array[j][1:m.K] .= swa.trKₖYstar_array[j] .- 2 .* swa.trKₖPYstar_array[j] .+ swa.trPKₖPYstar_array[j]
end

# function compute_trHHtD(m::cvc{TT, T}) where {TT, T}
#     dot(m.H, m.H̃)
# end

function update_trYstarV!(m::cvc{TT, T}) where {TT, T}
    # compute trHHᵀD where D = diag y_2^* - diag y_1^* y_1^*ᵀ
    trHHtD = dot(m.H, m.H̃)
    
    # compute y₁⋆ᵀSᵀSy₁⋆
    ystar1tHHtystar1 = dot(m.ystar1, m.HHtystar1)

    # compute trY⋆V where V = I - W(WᵀW)⁻¹Wᵀ
    sum!(m.trYstarV, m.ystar2)
    m.trYstarV .= m.trYstarV .- trHHtD .- ystar1tHHtystar1
end

function solve!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    @assert isassigned(m.trYstarV, 1)
    update_T!(m, swa, j)
    update_b!(m, swa, j)
    update_c!(m, swa, j)
    _solve!(m, swa, j)
end

function _solve!(m::cvc{TT, T}, swa::solve_working_array, j::Int) where {TT, T}
    swa.rvec_array[j][m.K + 1] = m.trYstarV[1]
    @views swa.gramMat_array[j][m.K + 1, 1:m.K] .= swa.gramMat_array[j][1:m.K, m.K + 1]
    @views swa.gramMat_array[j] .= Symmetric(transpose(swa.gramMat_array[j]))
    @views m.Jack[1:m.K+1, j] .= swa.gramMat_array[j] \ swa.rvec_array[j] # regular lsq
    @views m.Jack_nn[1:m.K+1, j] .= nonneg_lsq(swa.gramMat_array[j], swa.rvec_array[j], alg=:nnls, gram=true)
end

