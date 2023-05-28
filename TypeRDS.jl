############
### RDS Type
############

using Distributions, Intervals, Plots

""" The R

"""

struct RDS
    M::Interval                 # Phase Space M.
    SampleSpaceDimension::Int   # Dimesnion of Ω₀    
    LawOfSamples::Distribution  # Distribution of Ω₀
end

"""
    f<sub>ω</sub>   (ω::Float64, X)

Computes Xₖ₊₁ = fω(Xₖ).

## Returns
- Xₖ₊₁ = fω(Xₖ).

Note that 'fω' is a function modulo 1 so that Im(fω) ∈ M .

"""
function fω(ω::Float64, X)
    newvec = Vector{Float64}()
    for x in X
        push!(newvec, mod(ω * x, 1))
    end
    return newvec

end

"""
    struct PhaseSpaceDomainException <: Exception

An exception type representing an error related to the phase space domain.

## Fields
- `message::String`: A description of the exception.

"""
struct PhaseSpaceDomainException <: Exception
    message::String
end

"""
    sampleTraj(f::RDS, n::Int64, x0)

Returns sample trajectory of length n given an initial vector of data and random dynamical system.

## Fields
- `f': Random Dynamical System.
- 'n': Length of trajectory.
- 'x0': Initial data vector.
"""
function sampleTraj(f::RDS, n::Int64, x0) 
    if n <= 0
        throw(DomainError(n, "The number of iterations, $n, must be nonnegative."))
    end

    for x in x0
        if !(x in f.M)
            return  throw(PhaseSpaceDomainException("Initial data vector, $x0, must be in phase space M, but $x ∉ M."))
        end
    end

    traj = Vector{}()
    omegas = rand(f.LawOfSamples, n)    # [ω₁, ω₂, ⋯, ωₙ] 
    curr = x0

    for ω in omegas
        curr = fω(ω, curr)  # e.g. curr = x1 = f\_ω₁(x0)
        push!(traj, curr)   # add Xₖ to traj
    end

    return traj
end

