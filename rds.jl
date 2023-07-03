"""
This file provides functionality related to random dynamical systems (RDS) and trajectory sampling.

Usage:
1. Include this file in your Julia script or interactive session using `include("RDS.jl")`.
2. Make sure you have the necessary packages installed.

Functions:
- `RDS`: A struct representing a random dynamical system.
- `sampleTraj`: Generate a sample trajectory of initial data.
- `timeseries`: Generate a time series from an RDS using a given function.
- `empiricalAverage`: Compute the empirical average of a trajectory or time series.
"""

using Distributions, Intervals, Plots, StaticArrays, Statistics

mutable struct Domain
    dim::Int
    modulo_coordinates::Vector{Int}
end

mutable struct RDS
    M::Interval                 # Phase Space M.
    SampleSpaceDimension::Int   # Dimension of Ω₀
    LawOfSamples::Distribution  # Distribution of Ω₀
    func::Function               # fω (generic function)
end

"""
    struct PhaseSpaceDomainException <: Exception

An exception type representing an error related to the phase space domain.

## Fields:
- `message::String`: A description of the exception.

"""
struct PhaseSpaceDomainException <: Exception
    message::String
end

"""
    sampleTraj(f::RDS, n::Int64, x0)

Returns sample trajectory of length n given an initial vector of data and random dynamical system.

## Fields
- `system::RDS`: Random Dynamical System.
- `x0`: Initial data vector.
- `n::Int64`: Length of trajectory.

Make sure our field `func` for our RDS is defined :

    (ω, x) -> ⋯ and not (x , ω) -> ⋯

This is to assure our function fω works properly.
"""
function sampleTraj(system::RDS, x0, n::Int64; type="quenched", RO=false) 
    if n <= 0
        throw(DomainError(n, "The number of iterations, $n, must be nonnegative."))
    end

    for x in x0
        if !(x in system.M)
            return  throw(PhaseSpaceDomainException("Initial data vector, $x0, must be in phase space M, but $x ∉ M."))
        end
    end

    traj = Vector{}()
    curr = x0
    if type == "quenched"
        omegas = rand(system.LawOfSamples, n)    # [ω₁, ω₂, ⋯, ωₙ] 
        for ω in omegas
            curr = [mod(system.func(ω, x), 1) for x in curr]  # e.g. curr = x1 = f\_ω₁(x0)
            push!(traj, curr)   # add Xₖ to traj
        end
        if R0 == true
            return [SVector{n}(traj), SVector{n}(omegas)]
        end
    elseif type == "annealed"
        for i in 1:n
            newtraj = []
            curr = traj[i]
            omegas = rand(system.LawOfSamples, length(x0))
            for (j, x) in enumerate(curr)
                push!(newtraj, mod(system.func(omegas[j], x), 1))
            end
            push!(traj, newtraj)
        end
    end

    return SVector{n}(traj)
end

"""
    timeseries(traj::AbstractVector, ϕ::Function)

Compute a time series of data points using the given random dynamical system (`rds`), functions `fω`, `ϕ`, an initial state `x0`, and a number of iterations `n`.

## Arguments
- `traj::AbstractVector`: Trajectory.
- `ϕ::Function`: Function used to create our timeseries.


## Returns
- `timeseries`: A time series of transformed data points.

"""
function timeseries(traj::AbstractVector, ϕ::Function)
    return SVector{length(traj)}([ϕ.(pos) for pos in traj])
end

"""
    timeseries(traj::AbstractVector, omegas::AbstractVector, ϕ::Function)

Computes a Random Observable dependent on our vector of omegas, ϕ, and traj.

# Arguments
- `traj::AbstractVector`: A vector containing the trajectory data.
- `omegas::AbstractVector`: A vector of ω values.
- `ϕ::Function`: Function used to create our timeseries.

"""
function timeseries(traj::AbstractVector, omegas::AbstractVector, ϕ::Function)
    timeseries = Vector{}()
    for (i, pos) in enumerate(traj)
        push!(timeseries, ϕ.(omegas[i], pos))
    end

    return SVector{length(traj)}(timeseries)
end

"""
    empiricalAverage(traj::AbstractVector)

Compute the empirical average of a given trajectory.

## Arguments
- `traj::AbstractVector`: A trajectory represented as a vector.
"""
function empiricalAverage(traj::AbstractVector)
    tmp = zeros(length(traj[1]))
    
    for i in eachindex(traj)
        for j in eachindex(traj[1])
            tmp[j] += traj[i][j]
        end
    end
    
    return SVector{length(tmp)}(tmp / length(traj))
end

"""
    sampling(n::Int, distribution::Distribution)

Generate n valid samples from a specified distribution on the interval [0,1].

## Arguments
- `n::Int`: The number of samples to generate.
- `distribution::Distribution`: The distribution from which to generate the samples.
= `precision`: Key parameter determining precision of sample.
"""
function sampling(n::Int, distribution::Distribution; precision=false) # Slow when distribution=Normal()
    # If precision is want, we use BigFloat.
    if precision==true
         # To randomly sample BigFloats from distribution, we use the inverse CDF function.
        uni=rand(BigFloat, n)           
        samples=quantile(distribution, uni)
    else
        samples = rand(distribution, n)
    end
    transformed_samples = (samples .- minimum(samples)) / (maximum(samples) - minimum(samples))  # Transform samples to the interval [0, 1]
    valid_samples = filter(x -> 0 ≤ x ≤ 1, transformed_samples)  # Filter out values outside [lower, upper]
    return valid_samples
end