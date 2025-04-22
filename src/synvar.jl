### Code adapted from https://github.com/Hua-Zhou/CensoredLinearModels.jl ###

"""
    synthetic(ỹ, δ, km, ce; method = :KSvR)

Compute the synthetic variable of potentially right-censored time `ỹ`. 
`method=:KSvR` uses the method by Koul-Susarla-van Ryzin (1981) 
<https://doi.org/10.1214/aos/1176345644>; `method=:Leurgans` uses the method 
by Leurgans (1987) <https://doi.org/10.2307/2336144>.
"""
function synthetic(
    ỹ      :: Real,
    δ      :: Bool, 
    km     :: KaplanMeier{TT, TS},
    ce     :: Vector;
    moment :: Symbol = :first,
    method :: Symbol = :KSvR
    ) where {TT, TS}
    # last value in km.times less than or equal to T
    idx = searchsortedlast(km.times, ỹ)
    if moment == :first
        if method == :KSvR
            # Koul-Susarla-vanRyzin synthetic variable
            Ŝ  = ỹ ≤ km.times[1] ? one(TS) : km.survival[idx]
            sv = δ ? (ỹ / Ŝ) : zero(TS)
        elseif method == :Leurgans
            # Leurgans synthetic variable
            if ỹ ≤ km.times[1]
                sv = ỹ
            else
                sv = ce[idx] + (1 - km.survival[idx]) *
                    (ỹ - km.times[idx]) / km.survival[idx] + ỹ
            end
        else
            error("un-recognized method $method")
        end
    elseif moment == :second
        if method == :KSvR
            # Koul-Susarla-vanRyzin synthetic variable
            Ŝ  = ỹ ≤ km.times[1] ? one(TS) : km.survival[idx]
            sv = δ ? (ỹ^2 / Ŝ) : zero(TS)
        elseif method == :Leurgans
            # Leurgans synthetic variable
            if ỹ ≤ km.times[1]
                sv = ỹ^2
            else
                sv = ỹ^2 + ce[idx] + (1 - km.survival[idx]) * 
                    (ỹ^2 - km.times[idx]^2) / km.survival[idx]
            end
        else
            error("un-recognized method $method")
        end
    else
        error("un-supported moment $moment")
    end
    sv
end

"""
    synthetic!(out, ỹ, δ, km, ce; method = :KSvR)

Overwrite `out` by the synthetic variables, based on the observed times `ỹ` and 
censoring indicators `δ`. `method=:KSvR` uses the synthetic responses by 
Koul, Susarla, and Van Ryzin (1981) <https://doi.org/10.1214/aos/1176345644>; 
`method=:Leurgans` uses the method by Leurgans (1987) <https://doi.org/10.2307/2336144>.
"""
function synthetic!(
    out    :: AbstractVector,
    ỹ      :: AbstractVector,
    δ      :: Union{Vector{Bool}, BitVector, SubArray{Bool}},
    km     :: KaplanMeier,
    ce     :: AbstractVector;
    moment :: Symbol = :first,
    method :: Symbol = :KSvR
    )
    @inbounds for i in eachindex(out)
        out[i] = synthetic(ỹ[i], δ[i], km, ce, moment = moment, method = method)
    end
    out
end

function synthetic!(
    ystar1      :: AbstractVector,
    ystar2      :: AbstractVector,
    ỹ      :: AbstractVector,
    δ      :: Union{Vector{Bool}, BitVector, SubArray{Bool}},
    kmc     :: KaplanMeier,
    cec     :: AbstractVector,
    cec2     :: AbstractVector;
    first_moment_synvar   :: Symbol = :Leurgans,
    second_moment_synvar   :: Symbol = :Leurgans
    )
    if first_moment_synvar == :Leurgans || first_moment_synvar == :KSvR
        synthetic!(ystar1, ỹ, δ, kmc, cec; moment = :first, method = first_moment_synvar)
    else
        throw(ArgumentError("unrecognized method $(first_moment_synvar)"))
    end

    if second_moment_synvar == :Leurgans || second_moment_synvar == :KSvR
        synthetic!(ystar2, ỹ, δ, kmc, cec2; moment = :second, method = second_moment_synvar)
    else
        throw(ArgumentError("unrecognized method $(second_moment_synvar)"))
    end
end