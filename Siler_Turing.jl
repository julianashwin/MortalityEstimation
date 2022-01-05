"""
Use MCMC from Turing.jl to estimate Siler model on mortality data
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
function gompertz(a,b,c, ages)

    mort = exp.(a .+ b .* ages) .+ c

    return mort
end


m_dist= m_data[1]
plot(gompertz(-11, 0.11, 0.001, ages),ylim = (-0.1,0.5))
plot(siler(-2.8560, 1.7643, -7.3629, 0.0652, 0.0003, ages),ylim = (-0.1,0.5))
plot!(gompertz(-7.4391, 0.0659, 0.0015, ages),ylim = (-0.1,0.5))
scatter!(m_dist,ylim = (-0.1,0.5))


# Mean of IG is b/(a-1)
plot(InverseGamma(2, 0.01), xlim=(0,1))
plot(InverseGamma(2, 0.1), xlim=(0,1))
plot(InverseGamma(2, 6), xlim=(0,100))
plot(truncated(Normal(0.1, 0.1),0.0,1.0))

@model function gompertz_static(m_dist)
    # Our prior beliefs
    c0 ~ Normal(0, 10)
    c1 ~ InverseGamma(2, 0.1)
    d ~ InverseGamma(2, 0.01)
    σ ~ InverseGamma(2, 0.001)
    # The number of observations.
    N = length(m_dist)
    for nn in 1:N
        # Find mean using the gompertz mortality function
        m_mean = exp(c0 + c1 * (nn-1)) + d
        # Draw from truncated normal dist
        m_dist[nn] ~ truncated(Normal(m_mean, σ), 0.0, 1.0)
    end
end

iterations = 10000
chain_hmc = sample(gompertz_static(m_dist), HMC(ϵ, τ), iterations)
display(chain_hmc)
chain_smc = sample(gompertz_static(m_dist), SMC(), iterations)
display(chain_smc)
chain_nuts = sample(gompertz_static(m_dist), NUTS(0.65), iterations,chn=2)
display(chain_nuts)

chain = chain_nuts
histogram(chain[:c0])
histogram(chain[:c1])
histogram(chain[:d])
histogram(chain[:σ])


# Declare our Turing model for Siler
@model function siler_static(m_dist)
    # Our prior beliefs
    b0 ~ Normal(0, 10)
    b1 ~ InverseGamma(2, 0.1)
    c0 ~ Normal(0, 10)
    c1 ~ InverseGamma(2, 0.1)
    d ~ InverseGamma(2, 0.001)
    σ ~ InverseGamma(2, 0.001)

    # The number of observations.
    N = length(m_dist)
    for nn in 1:N
        # Find mean using the siler mortality function
        m_mean = exp(b0 - b1* (nn-1)) + exp(c0 + c1 * (nn-1)) + d
        # Draw from truncated normal dist
        m_dist[nn] ~ truncated(Normal(m_mean, σ), 0.0, 1.0)
    end
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.

ϵ = 0.05
τ = 10

# Start sampling.
chain_hmc = sample(siler_static(m_dist), HMC(ϵ, τ), iterations)
display(chain_hmc)
chain_smc = sample(siler_static(m_dist), SMC(), iterations)
display(chain_smc)
chain_nuts = sample(siler_static(m_dist), NUTS(0.65), iterations,chn=2)
display(chain_nuts)


# Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
chain = chain_nuts
histogram(chain[:b0])
histogram(chain[:b1])
histogram(chain[:c0])
histogram(chain[:c1])
histogram(chain[:d])
histogram(chain[:σ])
