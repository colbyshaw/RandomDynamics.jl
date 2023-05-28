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

#= This function will determine the operation of the function that we use:
    - 'a' or 'A' = Addition (ω + x)
    - 's1' or 'S1' = Subtraction Type 1 (ω - x)
    - 's2' or 'S2' = Subtraction Type 2 (x - ω)
    - 'm' or 'M' = Multiplication (ω * x)
    - 'd1' or 'D1' = Division Type 1 (ω / x)
    - 'd2' or 'D2' = Division Type 2 (x / ω)
    - 'e1' or 'E1' = Exponential Type 1 (ω ^ x)
    - 'e2' or 'E2' = Exponential Type 2 (x ^ ω)
    - More can be added here
=#
function operation(o::String, ω::Float64, x::Float64, rds::RDS)
    # I wish Julia had switch statement implementations...

    if o == "a" || o == "A"
        return mod(ω + x, rds.Mub)

    elseif o == "s1" || o == "S1"
        return mod(ω - x, rds.Mub)

    elseif o == "s2" || o == "S2"
        return mod(x - ω, rds.Mub)

    elseif o == "m" || o == "m"
        return mod(ω * x, rds.Mub)

    elseif o == "d1" || o == "D1"
        return mod(ω / x, rds.Mub)

    elseif o == "d2" || o == "D2"
        return mod(ω \ x, rds.Mub)

    # We have to be careful with the exponentials because we could get complex results (with negative numbers)
    elseif o == "e1" || o == "E1"
        check = Complex(ω) ^ x
        if imag(check) != 0
            println("Complex result! Try with a different distribution.")
            return 0
        else
            return mod(real(check), rds.Mub)
        end

    elseif o == "e2" || o == "E2"
        check = Complex(x) ^ ω
        if imag(check) != 0
            println("Complex result! Try with a different distribution.")
            return 0
        else
            return mod(real(check), rds.Mub)
        end

    else
        println("Unknown Operation! See documentation for implemented operations, or create your own!")
        return 0

    end
end


# Function applied to a sample and an input
# Note that the sampling is done now in sampleTraj instead
function fω(ω::Float64, x, rds::RDS, o::String)
    # Define the vector of values we will return
    returningVector = Vector{Float64}(undef, length(x))

    # Display the given parameters
    #println("M = [$(rds.Mlb), $(rds.Mub)]")
    #println("Dimension of Ω₀ = $(rds.Omega0Dim)")

    # # Need to check the dimension and the values of the input matches the phase space M, as x could be a vector
    # if rds.Omega0Dim != length(x)
    #     return "The dimensions of the input are incorrect! x is dimension $(length(x)) while Ω₀ is dimension $(rds.Omega0Dim)"
    # end

    # Want to print the type of distribution that "lawOfSamples" is for the RDS
    #println("Distribution: $(rds.lawOfSamples)")

    #= Create the method 
        f: Ω₀ × M → M, (ω, x) ↦ f\_ω (x) 
    =#
    # For now, we will let fω (x) = ω × x
    index = 1
    for val in x
        returningVector[index] = operation(o, ω, val, rds)
        #mod(ω * val, rds.Mub)
        #println(ω * x, mod(ω * val, rds.Mub))
        index += 1
    end
    
    # Return the outputing vector
    return returningVector
    
end