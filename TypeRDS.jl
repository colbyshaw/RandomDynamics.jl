############
### RDS Type
############

using Distributions, Intervals, Plots, StaticArrays, Observables

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
function fω(ω::Float64, X, func::Function)
    newvec = Vector{Float64}()
    for x in X
        push!(newvec, mod(func(ω, x), 1))
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
- 'func': f_ω
"""
function sampleTraj(f::RDS, n::Int64, x0, func::Function) 
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
        curr = fω(ω, curr, func)  # e.g. curr = x1 = f\_ω₁(x0)
        push!(traj, curr)   # add Xₖ to traj
    end

    return SVector{n}(traj)

end

function timeseries(rds::RDS, fω::Function, ϕ::Function, x0, n::Int)
    timeseries = Vector{}()
    initTraj = sampleTraj(rds, n, x0, fω)
    print(initTraj)
    for pos in initTraj
        tmp = Vector{Float64}()
        for data in pos
            push!(tmp, ϕ(data))
        end
        push!(timeseries, tmp)
    end
           
    return SVector{n}(timeseries)
end

#function empiricalAverage(rds::RDS, fω::Function, ϕ::Function, x0, n::Int)
#    tmp = zeros(length(x0))
#    data = timeseries(rds, f, ϕ, x0, n)

#    for i in eachindex(data)
#        for j in eachindex(data[i])
#            tmp[j] = data[i][j]
#        end
#    end
    
#    return SVector{length(x0)}(tmp / n)
#end

function empiricalAverage(traj::AbstractVector)
    tmp = zeros(length(traj[1]))
    
    for i in eachindex(traj)
        for j in eachindex(traj[1])
            tmp[j] += traj[i][j]
        end
    end
    
    return SVector{length(tmp)}(tmp / n)
end


##############
###Basic Tests
##############

rds = RDS(Interval{Closed, Closed}(0,1), 1, Normal())
ϕ(x) = x + (1-x) * .5
f(ω, x) = ω * x
d = rand(truncated(Normal(), 0, 1), 1)
#histogram(d)
n = 10

traj = timeseries(rds, f, ϕ, d, n)
empiricalAverage(traj)

