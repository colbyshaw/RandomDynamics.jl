
include("rds.jl")
include("graphics.jl")

pre= 2^8
setprecision(pre)

x=sampling(100000, Laplace(), precision="true")
histogram(x)

uni = rand(1000)
for i in eachindex(uni)
    println("element $i: $(uni[i])")
    x=quantile(Normal(), uni[i])
end

# Testing which univariate distributions can we sample BigFloat values using our sampling() as in "rds.jl"

univariate_distributions = [
    Bernoulli(), Binomial(), Cauchy(), Frechet(1,2), Exponential(), Geometric(), Gumbel(0,1), 
    Laplace(), NegativeBinomial(), Pareto(2,1), Poisson(), Rayleigh(1), Uniform(), Weibull(1,2)
    ]

    # No Beta, Gamma, InverseGamma, Chi, Chisq, FDist, TDist. Buffer error on LogNormal, Normal, Levy

    for dist in univariate_distributions
    sampling(1000000, dist, precision="true")
    println("Sampling process for $dist is complete!")
end
