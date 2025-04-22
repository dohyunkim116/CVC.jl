### Code adapted from https://github.com/Hua-Zhou/CensoredLinearModels.jl ###

"""
    KaplanMeier

An immutable type containing survivor function estimates computed
using the Kaplan-Meier method.

The type has the following fields:  
* `times`   : Distinct event times
* `nevents` : Number of observed events at each time
* `ncensor` : Number of right censored events at each time
* `natrisk` : Size of the risk set at each time
* `survival`: Estimate of the survival probability at each time
* `stderr`  : Standard error of the log survivor function at each time

Use `KaplanMeier(times, status)` to compute the estimates and construct
this type.
"""
struct KaplanMeier{TT <: Real, TS <: Real}
    # ordered unique times from smallest to largest
    times    :: Vector{TT}
    # number of failures at tᵢ
    nevents  :: Vector{Int}
    # number of censored subjects in [tᵢ, tᵢ₊₁)
    ncensor  :: Vector{Int}
    # number of subjects at risk for failing instantaneously prior to tᵢ,
    # including failures at tᵢ
    natrisk  :: Vector{Int}
    # survival estimate (right-continuous)
    survival :: Vector{TS}
    # standard error of log(survival)
    stderr   :: Vector{TS}
end

function KaplanMeier(
    times  :: AbstractVector{TT}, 
    status :: Union{Vector{Bool}, BitVector, SubArray{Bool}};
    TS :: Type = Float64,
    fit :: Bool = true
    ) where TT <: Real
    n  = length(times)
    km = KaplanMeier{TT, TS}(
        sizehint!(TT[], n),
        sizehint!(Int[], n), 
        sizehint!(Int[], n),
        sizehint!(Int[], n), 
        sizehint!(TS[], n), 
        sizehint!(TS[], n)
        )
    if fit
        p  = sortperm(times)
        t  = times[p]
        s  = status[p]
        _fit!(km, t, s)
    end
    km
end

"""
    _fit!(km, tte, status)

Kaplan-Meier estimator, assuming `tte` is sorted and does not include baseline 
time.
"""
function _fit!(
    km     :: KaplanMeier{TT, TS},
    tte    :: AbstractVector{TT},
    status :: Union{Vector{Bool}, BitVector}
    ) where {TT <: Real, TS <: Real}
    empty!(km.times)
    empty!(km.nevents)
    empty!(km.ncensor)
    empty!(km.natrisk)
    empty!(km.survival)
    empty!(km.stderr)
    nobs = length(tte)
    dᵢ = 0        # Number of observed events at time t
    cᵢ = 0        # Number of censored events at time t
    nᵢ = nobs     # Number remaining at risk at time t
    es = one(TS)  # Estimator starting point
    gw = zero(TS) # Standard Error starting point
    t_prev = typemin(TT)
    @inbounds for i in 1:nobs
        t = tte[i]
        s = status[i]
        # aggregate over tied times
        if t == t_prev
            dᵢ += s
            cᵢ += !s
            continue
        end
        # new distinct time
        if t_prev ≠ typemin(TT)
            if dᵢ > 0
                es *= 1 - dᵢ / nᵢ
                gw += dᵢ / (nᵢ * (nᵢ - dᵢ))
                push!(km.times   , t_prev)
                push!(km.nevents , dᵢ)
                push!(km.ncensor , cᵢ)
                push!(km.natrisk , nᵢ)
                push!(km.survival, es)
                push!(km.stderr  , sqrt(gw))
            else # dᵢ == 0
                isempty(km.ncensor) || (km.ncensor[end] += cᵢ)
            end
        end
        nᵢ -= dᵢ + cᵢ
        dᵢ = Int(s)
        cᵢ = Int(!s)
        t_prev = t
    end
    # We need to do this one more time to capture the last time
    # since everything in the loop is lagged
    push!(km.times   , t_prev)
    push!(km.nevents , dᵢ)
    push!(km.ncensor , cᵢ)
    push!(km.natrisk , nᵢ)
    push!(km.survival, es)
    push!(km.stderr  , sqrt(gw))
    km
end

length(km::KaplanMeier) = length(km.times)

"""
    survival(km, t)

Return probability of `T>t` based on the KaplanMeier estimate `km`.
"""
function survival(km::KaplanMeier, t::Real)
    idx = searchsortedlast(km.times, t)
    idx == 0 ? one(eltype(km.survival)) : km.survival[idx]
end

"""
    meanvar(km::KaplanMeier)

Calculate the mean and variance from Kaplan-Meier estimate of a distribution.
"""
function meanvar(km::KaplanMeier)
    μ  = km.times[1] * (1 - km.survival[1])
    σ² = km.times[1] * μ
    @inbounds for j in 2:length(km.times)
        px  = km.times[j] * (km.survival[j - 1] - km.survival[j])
        μ  += px
        σ² += km.times[j] * px
    end
    # if last survival probability > 0, assume it's 0
    if km.survival[end] > 0
        px  = km.times[end] * km.survival[end]
        μ  += px
        σ² += km.times[end] * px
    end
    σ² -= abs2(μ)
    μ, σ²
end

"""
    intgrl!(out, km)

Overwrite `out[i]` by the integral `∫ I(t ≤ km.times[i]) G(t) / (1 - G(t)) dt`. 
This is the term adding to the observed value `zᵢ` in Leurgan's method.
"""
function intgrl!(out::AbstractVector, km::KaplanMeier)
    # out[j] = ∫ I(km.times[j] ≥ t) G(t) / (1 - G(t)) dt
    out[1] = 0
    @inbounds for j in 2:length(km)
        # start integration only when survival probability is < 1
        if km.survival[j - 1] == 1
            out[j] = 0
        else
            out[j] = out[j - 1] + (1 - km.survival[j - 1]) *
                (km.times[j] - km.times[j - 1]) / km.survival[j - 1]
        end
    end
    out
end

intgrl(km::KaplanMeier) = intgrl!(similar(km.survival), km)

"""
    intgrl2!(out, km)

Overwrite `out[i]` by the integral `∫ I(t ≤ km.times[i]) 2tG(t) / (1 - G(t)) dt`. 
This is the term adding to the observed value `zᵢ²` in Leurgan's method.
"""
function intgrl2!(out::AbstractVector, km::KaplanMeier)
    # out[j] = ∫ I(km.times[j] ≥ t) 2t G(t) / (1 - G(t)) dt
    out[1] = 0
    @inbounds for j in 2:length(km)
        # start integration only when survival probability is < 1
        if km.survival[j - 1] == 1
            out[j] = 0
        else
            out[j] = out[j - 1] + (1 - km.survival[j - 1]) * 
                    (km.times[j]^2 - km.times[j - 1]^2) / km.survival[j - 1]
        end
    end
    out
end

intgrl2(km::KaplanMeier) = intgrl2!(similar(km.survival), km)