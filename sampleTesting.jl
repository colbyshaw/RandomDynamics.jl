###########
### Testing
### Thank you Colby :)
###########

include("secondaryObjectives.jl")

"""
    testing(x0::AbstractVector, distribution::Distribution, iterations::Int64)

Plots the trajectory of an initial vector for a given distribution.

## Arguments
- `x0`: An initial vector of data points.
- `distribution`: A vector of possible distributions for a sample space Ω.
- `iterations`: Length of trajectory.
"""
function testing(x0::AbstractVector, distribution::Distribution, iterations::Int64, func::Function)
    rds = RDS(0, 1, 1, distribution)
    data = sampleTraj(rds, iterations, x0, func)
    dataTraj = Vector{Vector{Float64}}(undef, length(x0))

    plot(title=distribution, xlabel="Iterations", ylabel="Values", legend=false, titlefontsize=10) # Format plot object.

    # Instantiates a vector of vectors to store trajectory of datapoints.
    for i in eachindex(dataTraj)
        dataTraj[i] = Float64[]
    end
        
    # Fills vector of trajectories of datapoints.
    for dataPoint in data
        for (i, point) in enumerate(dataPoint)
            push!(dataTraj[i], point)
        end
    end

    for i in eachindex(dataTraj)
        plot!(dataTraj[i])
    end
    display(current())
end

# Runner test (could be put in a "sampleRunner.jl" file eventually)
iterations = 100
function func(a, b)
    return a * b
end
samplesNum = 10
distList = [Normal(), Cauchy(), Laplace(), Gamma(), InverseGamma(), Beta(), Erlang(), LogNormal()]
for dist in distList
    x0 = rand(dist, samplesNum)
    for i in eachindex(x0)
        x0[i] = mod(x0[i], 1)
    end
    testing(x0, dist, iterations, func)
end

function φ(x)
    return x + 3
end

x1 = rand(Normal(), samplesNum)

randDynm = RDS(0, 1, 1, Normal())

timeseries(randDynm, φ, x1, samplesNum, func)
empiricalAverage(randDynm, φ, x1, samplesNum, func) # Seems to be close to 0.03 no matter the distribution (for large number of samples)!
#testing(x0, dist, iterations, func)