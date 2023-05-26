using Distributions, StaticArrays, Plots

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
    println("Samples Values: $samples")

    # Where the Markov chain progress will be tracked
    markovProgression = Vector{}()
    #markovProgression = Vector{}()
    push!(markovProgression, x0)
    iteration = x0
    for sample in samples
        iteration = fω(sample, iteration, rds)
        push!(markovProgression, iteration)
    end
    # index = 2
    # for sample in samples
    #     if index == length(markovProgression) + 1
    #         break
    #     end
    #     #println("Sample $i = $(samples[i])")
    #     #println(sample)
    #     nextX = fω(sample, x0, rds)
    #     #println(nextX)
    #     markovProgression[index] = nextX[1]
    #     index += 1
    # end
    # Would we want to store anything as a Static Array? We have to generate everything before storing because 
    return markovProgression
end

# Testing the example that Colby used
# Normal, Cauchy, Laplace, Gamma, InverseGamma
rds = RDS(0, 1, 1, Normal())
data = sampleTraj(rds, 250, [0.1, 0.9, 0.2, 0.7])

w = Vector{Float64}()
x = Vector{Float64}()
y = Vector{Float64}()
z = Vector{Float64}()
for sample in data
    push!(w, sample[1])
    push!(x, sample[2])
    push!(y, sample[3])
    push!(z, sample[4])
end

#println(x)
#println(y)

plot(w)
plot!(x)
plot!(y)
plot!(z)
#title!("$(rds.lawOfSamples) Distribution, 1000 Samples")
title!("$(rds.lawOfSamples) Distribution, 250 Samples", titlefontsize=12)
xlabel!("Iterations")
ylabel!("Markov Chain")