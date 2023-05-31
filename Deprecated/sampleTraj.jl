using Distributions, StaticArrays

include("typeRDS.jl")

# This is our main sampling method for our package using type RDS
function sampleTraj(f, distr, n, x0)
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end
    return 1
end

function sampleTraj(rds::RDS, n, x0)
    # Check if number of iterations inputted is sufficient
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end

    # Where the Markov chain progress will be tracked
    markovProgression = Vector{Float64}(undef, n)
    markovProgression[1] = x0
    for i in 2:n
        nextX = f(rds, i - 1)
    end
    # Would we want to store anything as a Static Array? We have to generate everything before storing because 
    return markovProgression
end

# Testing the example that Colby used
# Normal, Cauchy, Laplace, Gamma, InverseGamma
# rds = RDS(0, 1, 1, Laplace())

# input = rand(rds.lawOfSamples, 10)
# data = sampleTraj(rds, 100, input, operation::Function)
# #data = sampleTraj(rds, 100, [0.1, 0.9, 0.2, 0.7], operation::Function)

# w = Vector{Float64}()
# x = Vector{Float64}()
# y = Vector{Float64}()
# z = Vector{Float64}()
# for sample in data
#     push!(w, sample[1])
#     push!(x, sample[2])
#     push!(y, sample[3])
#     push!(z, sample[4])
# end

# #println(x)
# #println(y)
# plot(w)
# plot!(x)
# plot!(y)
# plot!(z)
# #title!("$(rds.lawOfSamples) Distribution, 1000 Samples")
# title!("$(rds.lawOfSamples) Distribution, 100 Samples", titlefontsize=12)
# xlabel!("Iterations")
# ylabel!("Markov Chain")