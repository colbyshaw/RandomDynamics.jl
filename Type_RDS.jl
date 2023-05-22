############
### RDS Type
############

using Distributions, Intervals

# Once declaring a random dynamical system, will it need to be changed over time?
# Probably not, so I am using struct rather than a mutable struct.

struct RDS
    M::Interval
    SampleSpaceDimension::Int
    LawOfSamples::Distribution
end

# Evalution of f_ω.

function f_ω()

end