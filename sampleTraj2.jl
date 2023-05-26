using Distributions, StaticArrays

include("typeRDS.jl")

function sampleTraj(f, distr, n, x0)
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end
    return 1
end

# This is our main sampling method for our package using type RDS
function sampleTraj(rds::RDS, n::Int64, x0)
    # Check if number of iterations inputted is sufficient
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end

    # Sample all of our values according to the distribution of the ω's 
    samples = rand(rds.lawOfSamples, n)
    #println("Samples Values: $samples")

    # Where the Markov chain progress will be tracked
    markovProgression = Vector{Float64}(undef, n)
    markovProgression[1] = x0
    for i in range(1, length(x0))
        println(samples[i])
        nextX = fω(samples[i], x0[i], rds);
        println(nextX)
    end
    # Would we want to store anything as a Static Array? We have to generate everything before storing because 
    return markovProgression
end