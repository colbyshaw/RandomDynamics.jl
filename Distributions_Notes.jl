#########################
### Distributions package
#########################

using Random, Distributions, Plots

# Setting seed
Random.seed!(123) 

# Normal() No parameters implies standard normal distribution
p(x) = 1/sqrt(2π) * exp(-x^2/2) # generic function 

# We can check paramaters of a univariate distribution d with params()
params(Normal())

# 100000 samples from above distribution
# rand(distribution, sample size) creates a new array. 
# rand!() takes an existing array and fills it with samples from the given distribution. For e.g.

# z = zeros(5)
# rand!(Normal(), z)


x = rand(Normal(), 1000000) 

histogram(x, normalize=:pdf, label=:"Sampled")
plot!(p, label=:"Std. Normal Dist.", color=:red)

# Fits best distribution to sample.
fit(Normal, x)

# Evaluate pdf of Normal() at 0.
pdf.(Normal(), 0)

# Evaluate cdf of Normal()
cdf.(Normal(), 0)

y = rand(Binomial(10, .55), 100000) 
histogram(y)

###########################
### Truncated Distributions
###########################

# Normal Distrbution truncated to [0, inf)
trunk = truncated(Normal(), lower=0)
histogram(rand(trunk, 100000), label="Truncated Normal Distribution", normalize=:pdf)

##############################
### Multivariate Distributions
##############################

μ = [0.0, 0.0]
Σ = [1.0 0.0; 0.0 1.0] 

# 10 samples.
x = zeros(2, 1)
mvn = MvNormal(μ, Σ)
sample = rand!(mvn, x)

# Note an easy way to write a matrix given data is by reshaping it.
matrix_1 = reshape([2,1,9,4,3,5], (2, 3))
matrix_2 = reshape([1,2,3], (3,1))

# Note: To convolute two probability distributions d1 and d2, use convolve(d1::Distribution, d2::Distribution).
