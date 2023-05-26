using Distributions

# Type Random Dynamical System used for our package
mutable struct RDS
    #= 
    Do we want to maybe write this so that we can use the "Range" package? This might be better though as
    we are then able to have endpoints to sample from instead of just some "Range" object.
    =#
    Mlb::Float64        # Lower bound of phase space M
    Mub::Float64        # Upper bound of phase space M
    Omega0Dim::Int64    # Dimension of Ω₀
    # lawOfSamples        # Law of samples of the ω values
    lawOfSamples::Distribution # If we only want to pull from the types of distributions from Distributions.jl (maybe write a different function?)
end

# Function applied to a sample and an input
# Note that the sampling is done now in sampleTraj instead
function fω(ω::Float64, x, rds::RDS)
    # Define the vector of values we will return
    returningVector = Vector{Float64}(undef, length(x))

    # Display the given parameters
    println("M = [$(rds.Mlb), $(rds.Mub)]")
    println("Dimension of Ω₀ = $(rds.Omega0Dim)")

    # Need to check the dimension and the values of the input matches the phase space M, as x could be a vector
    if rds.Omega0Dim != length(x)
        return "The dimensions of the input are incorrect! x is dimension $(length(x)) while Ω₀ is dimension $(rds.Omega0Dim)"
    end

    # Want to print the type of distribution that "lawOfSamples" is for the RDS
    println("Distribution: $(rds.lawOfSamples)")

    #= Create the method 
        f: Ω₀ × M → M, (ω, x) ↦ f\_ω (x) 
    =#
    # For now, we will let fω (x) = ω × x
    index = 1
    for val in x
        returningVector[index] = mod(ω * val, rds.Mub)
        #println(ω * x, mod(ω * val, rds.Mub))
        index += 1
    end
    
    # Return the outputing vector
    return returningVector
    
end