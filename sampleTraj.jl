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
    if n < 0
        return "n is $n, but n must be a positive number of iterations!"
    end
    samples = f(rds, x0)
    return samples
end