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

mutable struct RDSDomain
    dim::Int
    modulo_coordinates::Vector{Bool} 
end

mutable struct RDS
    M::Interval                 # Phase Space M.
    SampleSpaceDimension::Int   # Dimension of Ω₀
    LawOfSamples::Distribution  # Distribution of Ω₀
    func::Function               # fω (generic function)
end

"""
    makeDistribution(M::RDSDomain, f::Distribution, truncation::Int64)

Computes a new density, g, on our domain M given an initial distribuition f.

## Fields
- `M::RDSDomain`: Domain we are dealing with.
- `f`: Initial multivariate distribuition.
- `truncation::Int64`: Controls how accurate we want our new density, g, to be.

"""
function makeDistribution(M::RDSDomain, f::Distribution, truncation::Int64)
    if M.dim != length(f)
        throw(ArgumentError("Dimension for domains of M and f must match!"))
    end

    function g(x::Vector)
        range = -truncation:truncation
        total = pdf(f, x)
        zero = [i for i in 1:length(M.modulo_coordinates) if M.modulo_coordinates[i] == false]
        modulos = [i for i in 1:length(M.modulo_coordinates) if M.modulo_coordinates[i] == true]
        translator = Vector{}(undef, M.dim)

        for elt in Iterators.product(fill(range, length(modulos))...)
            # Take into account the modulo of our domain M.
            vals = collect(elt)
            translator[zero] .= 0
            translator[modulos] = vals
            if translator != zeros(M.dim)
                total += pdf(f, x .+ translator)
            end
        end
        return total
    end
    return g 
end

"""
    makeDistributionCross(distList)

Constructs the PDF of the distribution supported on the cross product of domains from the given distributions list.

## Arguments
- `distList`: List of distributions that each are defined on a given domain.

## Returns
- New described PDF for some new distribution.
"""
function makeDistributionCross(distList, weights)
    # Return a function that applies corresponding functions to a given input vector
    function h(x::Vector)
        returnVec = Vector{}()
        for (i, dist) in enumerate(distList)
            push!(returnVec, pdf(dist, x[i]))
        end
        
        combinedVals = sum(returnVec .* weights)                # Weighted sum of the PDF values
        normalizationConstant = sum(weights)                    # Sum of the weights
        adjustment = combinedVals / normalizationConstant       # Normalize the PDF to ensure it integrates to 1

        return adjustment
    end

    return h
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
        if RO == true
            return [SVector{n}(traj), SVector{n}(omegas)]
        end
    elseif type == "annealed"
        for i in 1:n
            newtraj = []
            if i != 1
                curr = traj[i-1]
            end
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
function timeseries(traj::AbstractVector, ϕ::Function, omegas::AbstractVector)
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
- `precision`: Key parameter determining precision of sample.
- 'BFNorm': Key parameter used if want BigFloat samples from the Standard Normal Distribution.
"""
function sampling(n::Int, distribution::Distribution; precision=false) 
    # If precision is true, we want to use BigFloat.
    if precision==true
        # To randomly sample BigFloats from distribution, we use the inverse CDF function.
        # Does not always work, depending on the given univariate distribution
        if isa(distribution, Normal)
            precision=2^9
            dist = Normal(BigFloat(0, precision), BigFloat(1, precision))
            samples = rand(dist, n)
        else
            uni=rand(BigFloat, n)
            samples=quantile(distribution, uni)
        end
    else
        samples = rand(distribution, n)
    end
    transformed_samples = (samples .- minimum(samples)) / (maximum(samples) - minimum(samples))  # Transform samples to the interval [0, 1]
    valid_samples = filter(x -> 0 ≤ x ≤ 1, transformed_samples)  # Filter out values outside [lower, upper]
    return valid_samples
end

"""
    makeBigFloat(num::Float64, digits::Int64)

Constructs a BigFloat number that approximates a given Float64 number by adding numbers up to the given digits amount.

## Arguments
- `num`: Some Float64 number (only about 14 decimal place precision).
- 'digits': Number of digits to add to the original Float64 number.

## Returns
-  BigFloat number.
"""
function makeBigFloat(num::Float64, digits::Int64)
    unif = Uniform()
    strNum = string(num)

    # Create blank array
    newNumArray = [""]
    for x in range(1, length(strNum) - 1)
        push!(newNumArray,"")
    end

    # Replace each array with the character position of the given number
    for i in range(1, length(strNum))
        newNumArray[i] = string(strNum[i])
    end

    addDigits = Vector{}()
    # Uniformly sample a digit from [0, 1, ..., 9] to add to the end of the sample
    for i in 1:digits
        push!(addDigits, string(Int(floor(10 * rand(unif)))))
    end

    # Add newly generated digits to the original number
    for val in addDigits
        push!(newNumArray, val)
    end

    # Make everything into one string
    newString = ""
    for dig in newNumArray
        newString *= dig
    end

    # Convert to BigFloat
    return parse(BigFloat, newString)
end