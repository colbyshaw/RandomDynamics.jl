using Distributions, Plots

# Set the desired precision
precision = 2^9  # Adjust as needed
# Be sure to check the different values for below!

# Create a BigFloat distribution
# dist1 = Beta(BigFloat(2, precision))                                                        # Correct Density and Incorrect BigFloat Values
# dist2 = Gamma(BigFloat(7.5, precision), BigFloat(1, precision))                             # Correct Density and Incorrect BigFloat Values
# dist3 = InverseGamma(BigFloat(3, precision), BigFloat(0.5, precision))                      # Correct Density and Incorrect BigFloat Values
# dist4 = Chi(BigFloat(1, precision))                                                         # Correct Density and Incorrect BigFloat Values
dist5 = Normal(BigFloat(0, precision), BigFloat(1, precision))

# No Beta, Gamma, InverseGamma, Chi, Chisq, FDist, TDist. Buffer error on LogNormal, Normal, Levy
# LogNormal, Normal, and Levy all use the Err function. There may be an integral error in the code for Err function. 

# Generate the samples
# samples = rand(dist1, 1000000)
# samples = rand(dist2, 1000000)
# samples = rand(dist3, 1000000)
# samples = rand(dist4, 1000000)
samples = rand(dist5, 1000000)

x = rand(samples)

println(x, "\n", typeof(x))

# Check the type of the samples
for val in samples
    if typeof(val) != BigFloat
        print("NOT BigFloat")
        break
    end
    # println(val)
end

# xLB, xUB, binsNum = 0, 1, 1000
# xLB, xUB, binsNum = 0, 20, 1000
# xLB, xUB, binsNum = 0, 1, 10000
# xLB, xUB, binsNum = 0, 3, 1000
xLB, xUB, binsNum = -4, 4, 1000

histogram(samples, bins=binsNum, xlims=(xLB, xUB))
display(current())