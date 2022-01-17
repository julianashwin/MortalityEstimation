"""
Use MCMC from Turing.jl to estimate Siler model on mortality data
"""

if occursin("jashwin", pwd())
    cd("C://Users/jashwin/Documents/GitHub/MortalityEstimation/")
else
    cd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
end


# Set up the multiple chains
using Distributed
Nchains = 2
#addprocs(Nchains)
#rmprocs(5)
workers()


@everywhere begin

    # Import libraries.
    using Turing, StatsPlots, Random, Optim, StatsBase, LinearAlgebra
    using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings

    # Import the mortality data
    mort_df = CSV.read("data/clean/USA_life.csv", DataFrame, ntasks = 1)
    showtable(mort_df)

    # Convert this into a matrix of mortality rates over time
    chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]
    m_data = chunk(mort_df.mx, 101)
    ages = Int64.(0:maximum(mort_df.Age))
    years = unique(mort_df.Year)

    @assert length(m_data)==length(years) "number of years doesn't match length of m_data"
    @assert length(m_data[1])==length(ages) "number of ages doesn't match length of m_data[1]"

end

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





@everywhere begin
    # Declare our Turing model for Siler
    @model function siler_static(m_dist, ages)
        # The number of observations.
        N = length(m_dist)
        # Our prior beliefs
        b0 ~ Normal(0, 10)
        b1 ~ InverseGamma(2, 0.1)
        c0 ~ Normal(0, 10)
        c1 ~ InverseGamma(2, 0.1)
        d ~ InverseGamma(2, 0.001)
        σ ~ InverseGamma(2, 0.001)
        # Variance matrix
        Σ = σ.*I(N)
        # Find mean using the siler mortality function
        m_means = exp.(b0 .- b1.* ages) .+ exp.(c0 .+ c1 .* ages) .+ d
        # Draw from truncated normal dist
        m_dist ~ MvNormal(m_means, Σ)

    end
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
ϵ = 0.05
τ = 10
iterations = 5000
# Sample for 1933
@time chain_1933 = pmap(x->sample(siler_static(m_data[1], ages), NUTS(0.65), iterations), 1:Nchains)

sample(siler_static(m_data[1], ages), NUTS(0.65), iterations,chn=4)
display(chain_1933)
# Sample for 2019
@time chain_2019 = pmap(x->sample(siler_static(m_data[87], ages), NUTS(0.65), iterations), 1:Nchains)
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
@model function siler_static_T(m_data)
    T = length(m_data)

    # Our prior beliefs
    b0 ~ filldist(Normal(0, 10), T)
    b1 ~ filldist(InverseGamma(2, 0.1), T)
    c0 ~ filldist(Normal(0, 10), T)
    c1 ~ filldist(InverseGamma(2, 0.1), T)
    d ~ filldist(InverseGamma(2, 0.001), T)
    σ ~ filldist(InverseGamma(2, 0.001), T)

    # The number of observations.
    for tt in 1:T
        N = length(m_data[tt])
        for nn in 1:N
            # Find mean using the siler mortality function
            m_mean = exp(b0[tt] - b1[tt]* (nn-1)) + exp(c0[tt] + c1[tt] * (nn-1)) + d[tt]
            # Draw from truncated normal dist
            m_data[tt][nn] ~ truncated(Normal(m_mean, σ[tt]), 0.0, 1.0)
        end
    end
end


@model function siler_dyn(m_data, ages)
    T = length(m_data)
    N = length(ages)

    b0 = Vector(undef, T)
    b1 = Vector(undef, T)
    c0 = Vector(undef, T)
    c1 = Vector(undef, T)
    d = Vector(undef, T)
    σ = Vector(undef, T)
    # Our prior beliefs
    b0[1] ~ Normal(0, 10)
    b1[1] ~ InverseGamma(2, 0.1)
    c0[1] ~ Normal(0, 10)
    c1[1] ~ InverseGamma(2, 0.1)
    d[1] ~ InverseGamma(2, 0.001)
    σ[1] ~ InverseGamma(2, 0.001)

    σ_b0 ~ InverseGamma(2, 0.1)
    σ_b1 ~ InverseGamma(2, 0.1)
    σ_c0 ~ InverseGamma(2, 0.1)
    σ_c1 ~ InverseGamma(2, 0.1)
    σ_d ~ InverseGamma(2, 0.1)
    σ_σ ~ InverseGamma(2, 0.1)

    # First period
    Σ = σ[1].*I(N)
    # Find mean using the siler mortality function
    m_means = exp.(b0[1] .- b1[1].* ages) .+ exp.(c0[1] .+ c1[1] .* ages) .+ d[1]
    # Draw from truncated normal dist
    m_data[1] ~ MvNormal(m_means, Σ)

    #filldist(Normal(0, 10), T)
    # The number of observations.
    for tt in 2:T
        # Update parameters
        b0[tt] ~  Normal(b0[tt-1], σ_b0)
        b1[tt] ~ truncated(Normal(b1[tt-1], σ_b1), 0.0, 1e10)
        c0[tt] ~ Normal(c0[tt-1], σ_c0)
        c1[tt] ~ truncated(Normal(c1[tt-1], σ_c1), 0.0, 1e10)
        d[tt] ~ truncated(Normal(d[tt-1], σ_d), 0.0, 1.0)
        σ[tt] ~ truncated(Normal(σ[tt-1], σ_σ), 1e-8, 1e10)
        # Variance matrix
        Σ = σ[tt].*I(N)
        # Find mean using the siler mortality function
        m_means = exp.(b0[tt] .- b1[tt].* ages) .+ exp.(c0[tt] .+ c1[tt] .* ages) .+ d[tt]
        # Draw from truncated normal dist
        m_data[tt] ~ MvNormal(m_means, Σ)
    end
end


periods = [1,10,20,30,40,50,60,70,80]
periods = [1,5]
years_selected = years[periods]
iterations = 100
chain_dyn = sample(siler_dyn(m_data[periods], ages), NUTS(0.65), iterations,chn=4)
display(chain_dyn)


b0s = Symbol.("b0[".*string.(1:length(periods)).*"]")
b1s = Symbol.("b1[".*string.(1:length(periods)).*"]")
c0s = Symbol.("c0[".*string.(1:length(periods)).*"]")
c1s = Symbol.("c1[".*string.(1:length(periods)).*"]")
ds = Symbol.("d[".*string.(1:length(periods)).*"]")
σs = Symbol.("σ[".*string.(1:length(periods)).*"]")

"""
Function to extract summary statistics from time series variables
"""
function summarise_stats(chain_in, years_selected)
    df_post = DataFrame(chain_in)
    rename!(df_post, vcat(names(df_post)[1:2], string.(years_selected)))
    df_ests = describe(df_post[:, 3:end])
    df_ests.variable = parse.(Int64,string.(df_ests.variable))
    df_ests[:,:pc975] .= 0.0
    df_ests[:,:pc025] .= 0.0
    for yy in 3:ncol(df_post)
        df_ests.pc975[(yy-2)] = percentile(df_post[:,yy], 97.5)
        df_ests.pc025[(yy-2)] = percentile(df_post[:,yy], 2.5)
    end

    return(df_ests)
end

# Extract the posterior distributions and some summary statistics
b0_ests = summarise_stats(chain_dyn[b0s], years_selected)
b1_ests = summarise_stats(chain_dyn[b1s], years_selected)
c0_ests = summarise_stats(chain_dyn[c0s], years_selected)
c1_ests = summarise_stats(chain_dyn[c1s], years_selected)
d_ests = summarise_stats(chain_dyn[ds], years_selected)
σ_ests = summarise_stats(chain_dyn[σs], years_selected)


plot(layout = (2,3))
plot!(b0_ests.variable, b0_ests.mean, title = L"b_{0}", label = false, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc975, linestyle = :dash, label = false, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc025, linestyle = :dash, label = false, subplot = 1)

plot!(b1_ests.variable, b1_ests.mean, title = L"b_{1}", label = false, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc975, linestyle = :dash, label = false, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc025, linestyle = :dash, label = false, subplot = 2)

plot!(c0_ests.variable, c0_ests.mean, title = L"c_{0}", label = false, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc975, linestyle = :dash, label = false, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc025, linestyle = :dash, label = false, subplot = 4)

plot!(c1_ests.variable, c1_ests.mean, title = L"c_{1}", label = false, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc975, linestyle = :dash, label = false, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc025, linestyle = :dash, label = false, subplot = 5)

plot!(d_ests.variable, d_ests.mean, title = L"d", label = false, subplot = 3)
plot!(d_ests.variable, d_ests.pc975, linestyle = :dash, label = false, subplot = 3)
plot!(d_ests.variable, d_ests.pc025, linestyle = :dash, label = false, subplot = 3)

plot!(σ_ests.variable, σ_ests.mean, title = L"\sigma", label = false, subplot = 6)#, ylim = (0.0, 0.1))
plot!(σ_ests.variable, σ_ests.pc975, linestyle = :dash, label = false, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc025, linestyle = :dash, label = false, subplot = 6)
savefig("figures/Siler_static/siler_dyn.pdf")







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
