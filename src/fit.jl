"""
    fit_me!(m::cvc)

# Positional arguments
- `m::cvc`: A `cvc` object.

# Keyword arguments
- `sla::Bool`: whether to use SnpLinAlg object for linear algera or not (default: false)
- `temp_dir_rm::Bool`: whether to remove the temporary directory or not (default: true)

# Output
- A concatenated vector of the following:
- `ϕ::Vector`: estimates of K + 1 components where K is the number of genetic components
- `h²::Vector`: a narrow-sense heritability estimate (one-element vector)
- `h²se::Vector`: block jacknife standard error estimate of h²
- `ϕ_nn::Vector`: estimates of K + 1 components with nonnegative constraints
- `h²_nn::Vector`: a narrow-sense heritability estimate with nonnegative constraints (one-element vector)
- `h²se_nn::Vector`: block jacknife standard error estimate of h² with nonnegative constraints
"""
function fit_me!(
    m          :: cvc{TT, T};
    sla        ::Bool = false,
    temp_dir_rm::Bool = true
    ) where {TT, T}
    for k in 1:m.K
        wa = working_array(m, k)
        for j in 1:m.J
            update_Xₖʲ!(m, wa, k, j, sla)
            update_L̃!(m, wa, k, j, sla)
            update_R̃!(m, wa, k, j, sla)
            update_S̃!(m, wa, k, j)
            update_F̃!(m, wa, k, j, sla)
            update_Ṽ!(m, wa, k, j, sla)
            update_U!(m, wa, k, j, sla)
            update_Ũ!(m, wa, k, j, sla)
            update_Q̃!(m, wa, k, j, sla)
        end
    end
    update_HHtU_k0!(m)

    # finalize the updates
    m.L̃ .+= m.L̃_0
    m.R̃ .+= m.R̃_0
    m.S̃ .+= m.S̃_0
    m.F̃ .+= m.F̃_0
    m.Ṽ .+= m.Ṽ_0

    for j in 1:m.J
        m.Q̃_kj_array[j] .+= m.Q̃_0
        m.U_kj_array[j] .+= m.U_0
        m.Ũ_kj_array[j] .+= m.Ũ_0
    end
    
    update_trYstarV!(m)
    swa = solve_working_array(m)
    # solve normal equations
    for j in 1:(m.J + 1)
        solve!(m, swa, j)
    end

    set_estimates!(m)
    m.isfitted[1] = true
    if temp_dir_rm
        rm.(readdir(m.temp_dir, join=true), recursive = true)
        #run(`rmdir $(m.temp_dir)`)
        #rm(m.temp_dir; recursive = true)
    end
    [m.ϕ; m.h2; m.h2se; m.ϕ_nn; m.h2_nn; m.h2se_nn]
end