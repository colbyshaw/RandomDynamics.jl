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
        #println(dataTraj[i])
    end
    display(current())
end

# Need to ensure that we can track each portion of the trajectory so that the histograms can all be tracked at each step (for graphing purposes)
# Do this by accessing the nth row and plotting it vs the actual PDF of the distribution the values were sampled from

# Want to create a testing function to show the histogram of each trajectory in sampleTraj
function testing2(x0::AbstractVector, distribution::Distribution, iterations::Int64, func::Function)
    rds = RDS(0, 1, 1, distribution)
    data = sampleTraj(rds, iterations, x0, func)
    dataTraj = Vector{Vector{Float64}}(undef, length(x0))

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

    # Graphs each of the histograms representing the values attained for the evolution of each starting value
    #plt = histogram(vals, label=:"Samples", normalize=:pdf, color=:blue) # Plot the values sampled above
    for i in eachindex(dataTraj)
        plt = histogram(dataTraj[i], label=:"Samples", normalize=:pdf, bins=0:0.01:1) # Plot the values sampled above
        display(plt)
    end
    #display(current())
end

# Runner test (could be put in a "sampleRunner.jl" file eventually)
iterations = 300
function func(a, b)
    return a * b
end
samplesNum = 300
#distList = [Normal(), Cauchy(), Laplace(), Gamma(), InverseGamma(), Beta(), Erlang(), LogNormal()]
distList = [Normal()]
for dist in distList
    x0 = rand(dist, samplesNum)
    for i in eachindex(x0)
        x0[i] = mod(x0[i], 1)
    end
    testing(x0, dist, iterations, func)
    testing2(x0, dist, iterations, func)
    
    vals = rand(dist, samplesNum)
    # for j in eachindex(vals)
    #     vals[j] = mod(vals[j], 1)
    # end
    plt = histogram(vals, label=:"Check", normalize=:pdf, bins=0:0.01:1) # Plot the values sampled above
    display(plt)
end

function φ(x)
    return x + 3
end

x1 = rand(Normal(), samplesNum)

randDynm = RDS(0, 1, 1, Normal())

timeseries(randDynm, φ, x1, samplesNum, func)
empiricalAverage(randDynm, φ, x1, samplesNum, func) # Seems to be close to 0.03 no matter the distribution (for large number of samples)!
#testing(x0, dist, iterations, func)