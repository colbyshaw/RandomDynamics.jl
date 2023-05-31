#= This function is an option for determining the operation of the function that we use:
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