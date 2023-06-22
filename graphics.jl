include("rds.jl")

using Statistics

"""
    testing(x0::AbstractVector, distribution::Distribution, iterations::Int64)

Plots the trajectory of an initial vector for a given distribution.

## Arguments
- `traj`: Trajectory we are testing.
"""

function testing(traj::AbstractVector)
    dataTraj = Vector{Vector{Float64}}(undef, length(traj[1]))

    plot(xlabel="Iterations", ylabel="Values", titlefontsize=10) # Format plot object.

    # Instantiates a vector of vectors to store trajectory of datapoints.
    for i in eachindex(dataTraj)
        dataTraj[i] = Float64[]
    end
        
    # Fills vector of trajectories of datapoints.
    for dataPoint in traj
        for (i, point) in enumerate(dataPoint)
            push!(dataTraj[i], point)
        end
    end

    for i in eachindex(dataTraj)
        plot!(dataTraj[i], legend=false)
    end
    display(current())
end

"""
    testing(x0::AbstractVector, distributions::Vector{Distribution{Univariate, Continuous}}, iterations::Int64)

Plots the trajectory of an initial vector of data for each distribution provided in 'distributions'.

## Arguments
- `x0`: An initial vector of data points.
- `distributions`: A vector of possible distributions for a sample space Ω.
- `func`: f_ω
- `iterations`: Length of trajectory.
"""
function testing(x0::AbstractVector, distributions::Vector{Distribution{Univariate, Continuous}}, func::Function, iterations::Int64)
    for dist in distributions
        rds = RDS(Interval{Closed, Closed}(0, 1), 1, dist)
        data = sampleTraj(rds, iterations, x0, func)
        dataTraj = Vector{Vector{Float64}}(undef, length(x0))

        plot(title=dist, xlabel="Iterations", ylabel="Values") # Format plot object.
        

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
            plot!(dataTraj[i], legend=false)
        end
        display(current())  # displays the plot associated with the distribution.
    end
end

"""
    tracking(traj::AbstractVector)

Tracks the evolution of the distributions of our states over time using a specified function, func (or fω).

## Argument
- `trajectory`: Trajectory we are displaying.

## Details
It is best for x0 to be sampled from a distribution.

The distribution evolution over time is recorded and displayed by generating a histogram for each iteration.
"""
function tracking(traj::AbstractVector)
    # Record and display distribution evolution over time.
    @gif for i in eachindex(traj)
        histogram(traj[i], xlims=(0,1), ylims=(0,10000), bins=50, xlabel="Values", ylabel="Frequency") # Plot the values sampled above
    end fps = 8E
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