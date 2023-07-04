"""
This file provides functionality related to random dynamical systems (RDS) and trajectory sampling.

Usage:
1. Include this file in your Julia script or interactive session using `include("TypeRDS.jl")`.
2. Make sure you have the necessary packages installed.

Functions:
- `RDS`: A struct representing a random dynamical system.
- `sampleTraj`: Generate a sample trajectory of initial data.
- `timeseries`: Generate a time series from an RDS using a given function.
- `empiricalAverage`: Compute the empirical average of a trajectory or time series.
"""

using Distributions, Intervals, Plots, StaticArrays

struct RDSDomain
    dim::Int
    modulo_coordinates::Vector{Int}
end

"""

Make sure 'func' is defined:

    f(ω, x) and not f(x, ω). 

This is to assure our function fω works properly.

"""
struct RDS
    M::Interval                 # Phase Space M.
    SampleSpaceDimension::Int   # Dimension of Ω₀    
    LawOfSamples::Distribution  # Distribution of Ω₀
    func::Function              # fω (generic function)
end

# Not using 'ObservedRDS' type anymore because of updated sampleTraj implementation

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
    sampleTraj(system::RDS, n::Int64, x0; type="quenched", RO=false)

Returns sample trajectory of length n given an initial vector of data and random dynamical system.

## Fields
- `system`: Random Dynamical System.
- `n`: Length of trajectory.
- `x0`: Initial data vector.
- 'type': Determines if dynamics will evolve according to either a quenched or annealed framework.
- 'RO': Determines if the list of ω values is returned to be used for Random Observables (RO)

"""
function sampleTraj(system::RDS, n::Int64, x0; type="quenched", RO=false)
    if n <= 0
        throw(DomainError(n, "The number of iterations, $n, must be nonnegative."))
    end

    for x in x0
        if !(x in system.M)
            return  throw(PhaseSpaceDomainException("Initial data vector, $x0, must be in phase space M, but $x ∉ M."))
        end
    end

    traj = Vector{}()
    push!(traj, x0)
    curr = x0

    if type == "quenched"
        omegas = rand(system.LawOfSamples, n)    # [ω₁, ω₂, ⋯, ωₙ]
        for ω in omegas
            curr = [mod(system.func(ω, x), 1) for x in curr]  # e.g. curr = x1 = f\_ω₁(x0)
            push!(traj, curr)   # add Xₖ to traj
        end
        if RO == true
            return [SVector{n+1}(traj), SVector{n}(omegas)]
        end
    elseif type == "annealed"
        for i in 1:n
            newtraj = []
            curr = traj[i]
            omegas = rand(system.LawOfSamples, length(x0))
            for (j, x) in enumerate(curr)
                push!(newtraj, mod(func(omegas[j], x), 1))
            end
            push!(traj, newtraj)
        end
    end

    return SVector{n+1}(traj)
end

"""
timeseries(traj::AbstractVector, ϕ::Function)

Compute a time series of data points using the given trajectory 'traj' and function 'ϕ' to evolve with.

## Arguments
- `traj`: Trajectory.
- `ϕ`: A function to apply to each data point.


## Returns
- `timeseries`: A time series of transformed data points.

"""
function timeseries(traj::AbstractVector, ϕ::Function)
    timeseries = Vector{}()
    count = 0               # Ensures the first initial state x0 does not have ϕ applied to it
    for pos in traj         # Change to enumerate as in function below
        if count == 0
            push!(timeseries, pos)
            continue
        end
        tmp = Vector{}()  # tmp = Vector{Float64}()
        for data in pos
            push!(tmp, ϕ(data))
        end
        # push!(timeseries, ϕ.(pos))   # See if ϕ can be applied to every step
        push!(timeseries, tmp)
        count += 1
    end
           
    return SVector{length(traj)}(timeseries)
end

"""
    timeseries(traj::AbstractVector, omegas::AbstractVector, ϕω::Function)

Compute a time series of data points using the given trajectory 'traj', omega values 'omegas' used, and the random observable function 'ϕω'.

## Arguments
- `traj`: Trajectory.
- 'omegas': Omega values used during the creation of the trajectory.
- `ϕω`: A function to apply to each data point using the ω values provided.


## Returns
- `timeseries`: A time series of transformed data points.

"""
function timeseries(traj::AbstractVector, omegas::AbstractVector, ϕω::Function)
    timeseries = Vector{}()
    for (i, pos) in enumerate(traj)
        if i == 1           # Ensures the first initial state x0 does not have ϕ applied to it
            push!(timeseries, pos)
            continue
        end
        tmp = Vector{}()  # tmp = Vector{Float64}()
        for data in pos
            push!(tmp, ϕω(omegas[i], data))
        end
        # push!(timeseries, ϕω.(pos))   # See if ϕ can be applied to every step
        push!(timeseries, tmp)
    end
           
    return SVector{length(traj)}(timeseries)
end

"""
    empiricalAverage(traj::AbstractVector)

Compute the empirical average of a given trajectory 'traj'.

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

"""
function sampling(n::Int, distribution::Distribution) # Make it capable to sample BigFloat
    samples = rand(distribution, n)  # Generate n samples from the specified distribution
    transformed_samples = (samples .- minimum(samples)) / (maximum(samples) - minimum(samples))  # Transform samples to the interval [0, 1]
    valid_samples = filter(x -> 0 ≤ x ≤ 1, transformed_samples)  # Filter out values outside [lower, upper]
    return valid_samples
end

# Testing for function


# function func(a, b)
#     return a * b
# end

# M = Interval{Closed, Closed}(0, 1)

# law = Normal()

# rds = RDS(M, 1, law, func)

# rds.func(-3, 12)