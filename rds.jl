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

struct RDS
    M::Interval                 # Phase Space M.
    SampleSpaceDimension::Int   # Dimension of Ω₀    
    LawOfSamples::Distribution  # Distribution of Ω₀
end

struct ObservedRDS <: RDS
    ϕω::Function                # Random Observable Function
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
    sampleTraj(system::RDS, n::Int64, x0, func::Function; type="quenched")

Returns sample trajectory of length n given an initial vector of data and random dynamical system.

## Fields
- `system`: Random Dynamical System.
- `n`: Length of trajectory.
- `x0`: Initial data vector.
- `func`: f_ω
- 'type': Determines if dynamics will evolve according to either a quenched or annealed framework.

Make sure 'func' is defined:

    f(ω, x) and not f(x, ω). 

This is to assure our function fω works properly.
"""
function sampleTraj(system::RDS, n::Int64, x0, func::Function; type="quenched") 
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
            curr = [mod(func(ω, x), 1) for x in curr]  # e.g. curr = x1 = f\_ω₁(x0)
            push!(traj, curr)   # add Xₖ to traj
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
    timeseries(rds::RDS, fω::Function, ϕ::Function, x0, n::Int)

Compute a time series of data points using the given random dynamical system (`rds`), functions `fω`, `ϕ`, an initial state `x0`, and a number of iterations `n`.

## Arguments
- `traj`: Trajectory.
- `ϕ`: A function to apply to each data point.


## Returns
- `timeseries`: A time series of transformed data points.

"""
function timeseries(traj::AbstractVector, ϕ::Function)
    timeseries = Vector{}()
    for pos in traj
        tmp = Vector{}()  # tmp = Vector{Float64}()
        for data in pos
            push!(tmp, ϕ(data))
        end
        push!(timeseries, tmp)
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
We now work to implement the timeseries and empiricalAverage functions above now with random observables ϕω.

For the timeseries, we must now be given a trajectory where the random noise was used along the trajectory. Thus, a new
function 'sampleTrajRO' (RO = Random Observables) is needed to do so.

This can be completed as below.
"""

"""
    sampleTrajRO(system::ObservedRDS, n::Int64, x0, func::Function; type="quenched")

Returns sample trajectory with respect to the random observables of length n given an initial vector of data and 
observed random dynamical system.

## Fields
- `system`: Observed Random Dynamical System.
- `n`: Length of trajectory.
- `x0`: Initial data vector.
- `func`: f_ω
- 'type': Determines if dynamics will evolve according to either a quenched or annealed framework.

Make sure 'func' is defined:

    f(ω, x) and not f(x, ω). 

This is to assure our function fω works properly.
"""
function sampleTrajRO(system::ObservedRDS, n::Int64, x0, func::Function; type="quenched") 
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
            curr = [system.ϕω(mod(func(ω, x), 1)) for x in curr]  # e.g. curr = X₁ = ϕω(f\_ω₁(x0))
            push!(traj, curr)   # add Xₖ to traj
        end
    elseif type == "annealed"
        for i in 1:n
            newtraj = []
            curr = traj[i]
            omegas = rand(system.LawOfSamples, length(x0))
            for (j, x) in enumerate(curr)
                push!(newtraj, system.ϕω(mod(func(omegas[j], x), 1)))   # Applies wrt ϕω
            end
            push!(traj, newtraj)
        end
    end

    return SVector{n+1}(traj)
end