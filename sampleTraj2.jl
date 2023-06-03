using Distributions, StaticArrays, Plots

include("typeRDS2.jl")

# For our current imple mentation, we would also need phase space (M) parameters
# Currently in progress (use the other function for now)
function sampleTraj(distr::Distribution, n::Int64, x0, func::Function, Mub::Float64, Mlb=0)
    # Check if number of iterations inputted is sufficient
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    elseif n == 0
        return x0
    end

    # Sample all of our values according to given distribution of the ω's
    samples = rand(distr, n)
    println("Samples Values: $samples")

    # Where the Markov chain progress will be tracked
    markovProgression = Vector{}()
    push!(markovProgression, x0)
    iteration = x0
    for sample in samples
        # Need to make an rds with the given distribution
        iteration = fω(sample, iteration, rds, func)
        push!(markovProgression, iteration)
    end
    # Would we want to store anything as a Static Array (SVector)?
    return markovProgression

end

# This is our main sampling method for our package using type RDS
function sampleTraj(rds::RDS, n::Int64, x0, func::Function)
    # Check if number of iterations inputted is sufficient
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    elseif n == 0
        return x0
    end

    # Sample all of our values according to the distribution of the ω's 
    samples = rand(rds.lawOfSamples, n)
    println("Samples Values: $samples")

    # Where the Markov chain progress will be tracked
    markovProgression = Vector{}()
    push!(markovProgression, x0)
    iteration = x0
    for sample in samples
        iteration = fω(sample, iteration, rds, func)
        push!(markovProgression, iteration)
    end
    # Would we want to store anything as a Static Array (SVector)?
    return markovProgression
end