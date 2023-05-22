############
### RDS Type
############

using Distributions, Intervals

# Once declaring a random dynamical system, will it need to be changed over time?
# Probably not, so I am using struct rather than a mutable struct.
struct RDS
    M = Interval{Closed, Closed}(0,1)
    SampleSpaceDimension::Int
    LawOfSamples::Distribution
end
