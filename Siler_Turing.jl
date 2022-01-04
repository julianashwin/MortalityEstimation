"""
Ue HMC from Turing.jl to estimate Siler model on mortality data
"""

if occursin("jashwin", pwd())
    cd("C://Users/jashwin/Documents/GitHub/MortalityEstimation/")
else
    cd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
end

# Import libraries.
using Turing, StatsPlots, Random
using CSV, DataFrames, TableView, StatsPlots

# Import the mortality data
mort_df = CSV.read("data/clean/USA_life.csv", DataFrame, ntasks = 1)
showtable(mort_df)

# Convert this into a matrix of mortality rates over time
chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]
m_data = chunk(m, 101)
ages = Int64.(0:maximum(mort_df.Age))
years = unique(mort_df.Year)

@assert length(m_data)==length(years) "number of years doesn't match length of m_data"
@assert length(m_data[1])==length(ages) "number of ages doesn't match length of m_data[1]"

# Set the true probability of heads in a coin.
function siler(b0,b1,c0,c1,d, ages)

    mort = exp.(b0 .- b1.* ages) .+ exp.(c0 .+ c1 .* ages) .+ d

    return mort
end

plot(siler(0.7494, 5.1240, 0.2916, 1.5719, 0.3386, ages))

data = rand(Bernoulli(p_true), last(Ns))

plot(InverseGamma(2, 3),xlim =(0,100))
plot(Normal(0, 0.01))


m_dist= m_data[1]
# Declare our Turing model.
@model function siler_static(m_dist)
    # Our prior belief about the probability of heads in a coin.
    b0 ~ InverseGamma(2, 3)
    b1 ~ InverseGamma(2, 3)
    c0 ~ InverseGamma(2, 3)
    c1 ~ InverseGamma(2, 3)
    d ~ InverseGamma(2, 3)

    # The number of observations.
    N = length(m_dist)
    for nn in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        m_dist[nn] ~ exp(b0 - b1* (nn-1)) + exp(c0 + c1 * (nn-1)) + d + Normal(0,0.01)
    end
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 10000
ϵ = 0.05
τ = 10

# Start sampling.
chain = sample(siler_static(m_dist), HMC(ϵ, τ), iterations)
chain = sample(siler_static(m_dist), SMC(), iterations)

# Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
histogram(chain[:b0])
histogram(chain[:b1])
histogram(chain[:c0])
histogram(chain[:c1])
histogram(chain[:d])
