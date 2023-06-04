###########
### Secondary Objectives
### Thank you Colby :)
###########

include("sampleTraj2.jl")

# Tracks the evolution of a time series with respect to the function φ
function timeseries(rds::RDS, φ::Function, x0, n::Int64, func::Function)
    timeseries = Vector{}()
    initTraj = sampleTraj(rds, n, x0, func) # Sample Trajectory

    for pos in initTraj
        tmp = Vector{Float64}()
        for data in pos
            push!(tmp, φ(data))
        end
        push!(timeseries, tmp)
    end

    return timeseries
end

# Averages the time series function values
function empiricalAverage(rds::RDS, φ::Function, x0, n::Int64, func::Function)
    tmp = zeros(size(x0))
    data = timeseries(rds, φ, x0, n, func)
    for i in 1:n
        for j in eachindex(data)
            tmp[i] = data[j][i]
        end
    end

    return SVector{length(x0)}(tmp / n)[1]
end