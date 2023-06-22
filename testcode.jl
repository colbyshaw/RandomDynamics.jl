
include("rds.jl")
include("graphics.jl")

using IntervalArithmetic

pre=2^7
setprecision(pre)

x0 = sampling(1000, Normal())
x0 = BigFloat.(x0)
#x0 = interval.(x0)

f(ω, x) = mod(ω, 1) * x

rds = RDS(Interval{Closed, Closed}(0,1), 1, Normal())
sampleTraj(rds, 10, x0, f)

