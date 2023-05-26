############
### RDS Type
############

using Distributions, Intervals, Plots

# Once declaring a random dynamical system, will it need to be changed over time?
# Probably not, so I am using struct rather than a mutable struct.

struct RDS
    M::Interval
    SampleSpaceDimension::Int
    LawOfSamples::Distribution
end

# Evalution of f_ω(x).
function fω(ω::Float64, X)
    newvec = Vector{Float64}()
    for x in X
        push!(newvec, mod(ω * x, 1))
    end
    return newvec

end

function sampleTraj(f::RDS, n::Int64, x0)
    if n <= 0
        return "Parameter n needs to be greater than 0."
    end

    for x in x0
        if !(x in f.M)
            return "x0 must be in our phase space M."
        end
    end

    traj = Vector{}()
    omegas = rand(f.LawOfSamples, n)
    curr = x0

    for ω in omegas
        curr = fω(ω, curr)
        push!(traj, curr)
    end

    return traj
end

rds = RDS(Interval{Closed, Closed}(0, 1), 1, Normal())
data = sampleTraj(rds, 20, [.9,.1])

x = Vector{Float64}()
y = Vector{Float64}()
for sample in data
    push!(x,sample[1])
    push!(y, sample[2])
end

plot(x)
plot!(y)