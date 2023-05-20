using Distributions

mutable struct RDS
    #= 
    Do we want to maybe write this so that we can use the "Range" package? This might be better though as
    we are then able to have endpoints to sample from instead of just some "Range" object.
    =#
    Mlb::Float64        # Lower bound of phase space M
    Mub::Float64        # Upper bound of phase space M
    Omega0Dim::Int64    # Dimension of Ω₀
    lawofSamples        # Law of samples of the ω values
    # lawofSamples::Distribution # If we only want to pull from the types of distributions from Distributions.jl
end

function MarkovGenerator(rds::RDS)
    println("M = [$(rds.Mlb), $(rds.Mub)]")
    println("Dimension of Ω₀ = $(rds.Omega0Dim)")
    # Want to print the type of distribution that "lawofSamples" is for the RDS
    #= Next, we want to create the method 
        f: Ω₀ × M → M, (ω, x) ↦ f\_ω (x) 
    =#
end