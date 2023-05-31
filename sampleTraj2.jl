using Distributions, StaticArrays, Plots

include("typeRDS2.jl")

function sampleTraj(f, distr, n, x0)
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end
    return 1
end

# This is our main sampling method for our package using type RDS
function sampleTraj(rds::RDS, n::Int64, x0, func::Function)
    # Check if number of iterations inputted is sufficient
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
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