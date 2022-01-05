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
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings

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


"""
Plot priors
"""
# Mean of IG is b/(a-1)
plot(layout = (2,3), yticks = false, size = (800,400))
plot!(Normal(0, 10), title = L"b_{0,t} \sim \mathcal{N}(0,10)",
    label = false, subplot = 1)
plot!(InverseGamma(2, 0.1), xlim=(0,1), title = L"b_{1,t} \sim \mathcal{IG}(2,0.1)",
    label = false, subplot = 2)
plot!(Normal(0, 10), title = L"c_{0,t} \sim \mathcal{N}(0,10)",
    label = false, subplot = 4)
plot!(InverseGamma(2, 0.1), xlim=(0,1), title = L"c_{1,t} \sim \mathcal{IG}(2,0.1)",
    label = false, subplot = 5)
plot!(InverseGamma(2, 0.001), title = L"d_{t} \sim \mathcal{IG}(2,0.001)",
    label = false, subplot = 3)
plot!(InverseGamma(2, 0.001), title = L"\sigma_{t} \sim \mathcal{IG}(2,0.001)",
    label = false, subplot = 6)
savefig("figures/Siler_static/siler_priors.pdf")





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

# Sample for 1933
chain_1933 = sample(siler_static(m_data[1]), NUTS(0.65), iterations,chn=4)
display(chain_1933)
# Sample for 2019
chain_2019 = sample(siler_static(m_data[87]), NUTS(0.65), iterations,chn=4)
display(chain_2019)

# Plot model fit
plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
scatter!(m_data[1], markershape = :cross, markeralpha = 0.5, label = "Data (1933)")
plot!(siler(median(chain_1933[:b0]), median(chain_1933[:b1]), median(chain_1933[:c0]),
    median(chain_1933[:c1]), median(chain_1933[:d]), ages), label = "Siler MCMC fit (1993)")
scatter!(m_data[87], markershape = :xcross, markeralpha = 0.5, label = "Data (2019)")
plot!(siler(median(chain_2019[:b0]), median(chain_2019[:b1]), median(chain_2019[:c0]),
    median(chain_2019[:c1]), median(chain_2019[:d]), ages), label = "Siler MCMC fit (2019)")
savefig("figures/Siler_static/siler_fit.pdf")



plot(siler(-2.8560, 1.7643, -7.3629, 0.0652, 0.0003, ages))

# Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
plot(layout = (2,3))
density!(chain_1933[:b0], title = L"b_{0}", label = "1933", subplot = 1)
density!(chain_2019[:b0], label = "2019", subplot = 1, legend = :topleft)
density!(chain_1933[:b1], title = L"b_{1}", subplot = 2, legend = false)
density!(chain_2019[:b1], subplot = 2, legend = false)
density!(chain_1933[:c0], title = L"c_{0}", subplot = 4, legend = false)
density!(chain_2019[:c0], subplot = 4, legend = false)
density!(chain_1933[:c1], title = L"c_{1}", subplot = 5, legend = false)
density!(chain_2019[:c1], subplot = 5, legend = false)
density!(chain_1933[:d], title = L"d", subplot = 3, legend = false)
density!(chain_2019[:d], subplot = 3, legend = false)
density!(chain_1933[:σ], title = L"\sigma", subplot = 6, legend = false)
density!(chain_2019[:σ], subplot = 6, legend = false)
savefig("figures/Siler_static/siler_1933vs2019.pdf")



# Declare our Turing model for dynamic Siler equation
@model function siler_dyn(m_data)
    T = length(m_data)

    # Our prior beliefs
    b0 ~ fillldist(Normal(0, 10))
    b1 ~ InverseGamma(2, 0.1)
    c0 ~ Normal(0, 10)
    c1 ~ InverseGamma(2, 0.1)
    d ~ InverseGamma(2, 0.001)
    σ ~ InverseGamma(2, 0.001)

    # The number of observations.
    for tt in 1:T
        m_dist = m_data[tt]
        N = length(m_dist)
        for nn in 1:N
            # Find mean using the siler mortality function
            m_mean = exp(b0 - b1* (nn-1)) + exp(c0 + c1 * (nn-1)) + d
            # Draw from truncated normal dist
            m_dist[nn] ~ truncated(Normal(m_mean, σ), 0.0, 1.0)
        end
    end
end




















##################################################################################

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
