include("rds.jl")

"""
    testing(x0::AbstractVector, distribution::Distribution, iterations::Int64)

Plots the trajectory of an initial vector for a given distribution.

## Arguments
- `traj::AbstractVector`: Trajectory we are testing.
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
    tracking(traj::AbstractVector)

Tracks the evolution of the distributions of our states over time using a specified function, func (or fω).

## Argument
- `trajectory::AbstractVector`: Trajectory we are displaying.

## Details
It is best for x0 to be sampled from a distribution.

The distribution evolution over time is recorded and displayed by generating a histogram for each iteration.
"""
function tracking(traj::AbstractVector)
    # Record and display distribution evolution over time.
    @gif for i in eachindex(traj) 
        histogram(traj[i], xlims=(0,1), bins=100, xlabel="Values", ylabel="Frequency", legend=false) 
    end fps = 8
end
   
function distDisplay(M::RDSDomain, f::Distribution, xrange::AbstractRange, yrange::AbstractRange, truncation::Int64)
    @gif for i in 1:truncation
        g = makeDistribution(M, f, i)
        zs = [g([xs,ys]) for ys in yrange, xs in xrange]
        surface(x, y, zs, xlabel="X", ylabel="Y", zlabel="Z", title="Density on M.")
    end fps = 2
end

dom = RDSDomain(2, [true, false])
f = MultivariateNormal([0, 0], [1 0.0; 0.0 1])
# Generate grid points for plotting
x = 0:.1:1
y = -4:.1:4

distDisplay(dom , f, x, y, 5)