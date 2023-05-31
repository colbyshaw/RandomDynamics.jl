###########
### Testing
### Thank you Colby :)
###########

include("typeRDS2.jl")
include("sampleTraj2.jl")

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
dist = Cauchy()
iterations = 100
function func(a, b)
    return a * b
end
samplesNum = 100
x0 = rand(dist, samplesNum)

testing(x0, dist, iterations, func)


"""
    testing(x0::AbstractVector, distributions::Vector{Distribution{Univariate, Continuous}}, iterations::Int64)

Plots the trajectory of an initial vector of data for each distribution provided in 'distributions'.

## Arguments
- `x0`: An initial vector of data points.
- `distributions`: A vector of possible distributions for a sample space Ω.
- `iterations`: Length of trajectory.
"""
function testing(x0::AbstractVector, distributions::Vector{Distribution{Univariate, Continuous}}, iterations::Int64)
    for dist in distributions
        rds = RDS(Interval{Closed, Closed}(0, 1), 1, dist)
        data = sampleTraj(rds, iterations, x0)
        dataTraj = Vector{Vector{Float64}}(undef, length(x0))

        plot(title=dist, xlabel="Iterations", ylabel="Values", titlefontsize=10) # Format plot object.
        

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
        display(current())  # displays the plot associated with the distribution.
    end
end
