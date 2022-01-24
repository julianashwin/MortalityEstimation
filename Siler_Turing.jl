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
#Nchains = 2
#addprocs(Nchains)
#rmprocs(5)
workers()


# Import libraries.
using Turing, StatsPlots, Random, Optim, StatsBase, LinearAlgebra, Optim
using TruncatedDistributions
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings

include("helpers.jl")


# Import the mortality data
mort_df = CSV.read("data/clean/SWE_life.csv", DataFrame, ntasks = 1)
showtable(mort_df)

# Convert this into a matrix of mortality rates over time
chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]
m_data = chunk(mort_df.mx, 101)
ages = Int64.(0:maximum(mort_df.Age))
years = unique(mort_df.Year)
T = length(m_data)

@assert length(m_data)==length(years) "number of years doesn't match length of m_data"
@assert length(m_data[1])==length(ages) "number of ages doesn't match length of m_data[1]"

# Set the true probability of heads in a coin.
function siler(B,b,C,c,d, ages)

    mort = exp.(- b.* (ages .+ B)) .+ exp.(c .* (ages.- C)) .+ d

    return mort
end


m_dist= m_data[1]
plot(siler(1.5, 1.7643, 113.0, 0.0652, 0.0003, ages),ylim = (-0.1,0.6))
scatter!(m_dist,ylim = (-0.1,0.6))
scatter!(m_data[T],ylim = (-0.1,0.6))


"""
Plot priors
"""
# Mean of IG is b/(a-1)
plot(layout = (2,3), yticks = false, size = (800,400))
plot!(InverseGamma(2, 10), title = L"B_{t} \sim \mathcal{IG}(2,10)",
    label = false, subplot = 1, xlims = (0,100))
plot!(InverseGamma(2, 1.0), xlim=(0,10), title = L"b_{t} \sim \mathcal{IG}(2,1)",
    label = false, subplot = 2)
plot!(InverseGamma(2, 80), title = L"C_{t} \sim \mathcal{IG}(2,10)",
    label = false, subplot = 4, xlims = (0,100))
plot!(InverseGamma(2, 0.1), xlim=(0,1), title = L"c_{t} \sim \mathcal{IG}(2,0.1)",
    label = false, subplot = 5)
plot!(InverseGamma(2, 0.0001), title = L"d_{t} \sim \mathcal{IG}(2,0.0001)",
    label = false, subplot = 3)
plot!(InverseGamma(2, 0.01), title = L"\sigma_{t} \sim \mathcal{IG}(2,0.01)",
    label = false, subplot = 6)
savefig("figures/Siler_static/siler_priors.pdf")


plot(LogNormal(log(0.01), 0.5))
mean(LogNormal(log(0.01), 0.5))

"""
Static Siler model
"""
# Declare our Turing model for Siler
@model function siler_static(m_dist, ages)
    # The number of observations.
    N = length(m_dist)
    # Our prior beliefs
    B ~ InverseGamma(2, 10)
    b ~ InverseGamma(2, 0.1)
    C ~ InverseGamma(2, 80)
    c ~ InverseGamma(2, 0.1)
    d ~ InverseGamma(2, 0.001)
    σ ~ InverseGamma(2, 0.001)
    # Variance matrix
    Σ = σ.*I(N)
    # Find mean using the siler mortality function
    m_means = exp.(- b.* (ages .+ B)) .+ exp.(c .* (ages .- C)) .+ d
    # Draw from truncated normal dist
    m_dist ~ MvNormal(m_means, Σ)

end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
ϵ = 0.05
τ = 10
iterations = 10000
# Sample for 1933
#@time chain_1933 = pmap(x->sample(siler_static(m_data[1], ages), NUTS(0.65), iterations), 1:Nchains)

map_static_1 = optimize(siler_static(m_data[1], ages), MAP())

@time chain_1 = sample(siler_static(m_data[1], ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_static_1.values.array)
display(chain_1)
# Sample for 2019
map_static_T = optimize(siler_static(m_data[T], ages), MAP())
@time chain_T = sample(siler_static(m_data[T], ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_static_T.values.array)
display(chain_T)

# Plot model fit
plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
scatter!(m_data[1], markershape = :cross, markeralpha = 0.5, label = "Data ("*string(years[1])*")")
plot!(siler(median(chain_1[:B]), median(chain_1[:b]), median(chain_1[:C]),
    median(chain_1[:c]), median(chain_1[:d]), ages), label = "Siler MCMC fit ("*string(years[1])*")")
scatter!(m_data[T], markershape = :xcross, markeralpha = 0.5, label = "Data ("*string(years[T])*")")
plot!(siler(median(chain_T[:B]), median(chain_T[:b]), median(chain_T[:C]),
    median(chain_T[:c]), median(chain_T[:d]), ages), label = "Siler MCMC fit ("*string(years[T])*")")
savefig("figures/Siler_static/siler_fit.pdf")

# Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
plot(layout = (2,3), size = (800, 400))
density!(chain_1933[:B], title = L"B", label = "1933", subplot = 1)
density!(chain_2019[:B], label = "2019", subplot = 1, legend = :topright)
density!(chain_1933[:b], title = L"b", subplot = 2, legend = false)
density!(chain_2019[:b], subplot = 2, legend = false)
density!(chain_1933[:C], title = L"C", subplot = 4, legend = false)
density!(chain_2019[:C], subplot = 4, legend = false)
density!(chain_1933[:c], title = L"c", subplot = 5, legend = false)
density!(chain_2019[:c], subplot = 5, legend = false)
density!(chain_1933[:d], title = L"d", subplot = 3, legend = false)
density!(chain_2019[:d], subplot = 3, legend = false)
density!(chain_1933[:σ], title = L"\sigma", subplot = 6, legend = false)
density!(chain_2019[:σ], subplot = 6, legend = false, xrotation = 45.0, margin=3Plots.mm)
savefig("figures/Siler_static/siler_1933vs2019.pdf")


"""
Multiple independent Siler models
"""
# Define Turing model
@model function siler_indep(m_data, ages)
    # Dimensions
    T = length(m_data)
    N = length(ages)
    # Our prior beliefs
    B ~ filldist(InverseGamma(2, 10), T)
    b ~ filldist(InverseGamma(2, 0.1), T)
    C ~ filldist(InverseGamma(2, 80), T)
    c ~ filldist(InverseGamma(2, 0.1), T)
    d ~ filldist(InverseGamma(2, 0.001), T)
    σ ~ filldist(InverseGamma(2, 0.001), T)

    # The number of observations.
    for tt in 1:T
        Σ = σ[tt].*I(N)
        # Find mean using the siler mortality function
        m_means = exp.(-b[tt].*(ages .+ B[tt])) .+ exp.(c[tt].*(ages .- C[tt])) .+ d[tt]
        # Draw from truncated normal dist
        m_data[tt] ~ MvNormal(m_means, Σ)
        # Find mean using the siler mortality function
    end
end

# Estimate the model
periods = Int.(1:20:T)
years_selected = years[periods]
iterations = 2000
# Find MAP estimate as starting point
map_indep = optimize(siler_indep(m_data[periods], ages), MAP(), AcceleratedGradientDescent(),
    Optim.Options(iterations=10_000, allow_f_increases=true))
# MCMC sampling
@time chain_indep = sample(siler_indep(m_data[periods], ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_indep.values.array)
display(chain_indep)

parests_indep = extract_variables(chain_indep, periods, years_selected)

plot_siler_params(parests_indep)
savefig("figures/Siler_dynamic/siler_staticT.pdf")






"""
Dynamic Siler model
"""
# Declare our Turing model for dynamic Siler equation
@model function siler_dyn(m_data, ages)
    # Dimensions
    T = length(m_data)
    N = length(ages)
    # Priors on variance terms
    σ_b0 ~ InverseGamma(2, 0.1)
    σ_b1 ~ InverseGamma(2, 0.01)
    σ_c0 ~ InverseGamma(2, 0.1)
    σ_c1 ~ InverseGamma(2, 0.01)
    σ_d ~ InverseGamma(2, 0.0001)
    σ_σ ~ InverseGamma(2, 0.001)
    # First period fixed
    b0 = Vector(undef, T)
    b1 = Vector(undef, T)
    c0 = Vector(undef, T)
    c1 = Vector(undef, T)
    d = Vector(undef, T)
    σ = Vector(undef, T)
    # First period from priors
    b0[1] ~ InverseGamma(2, 10)
    b1[1] ~ InverseGamma(2, 10)
    c0[1] ~ InverseGamma(2, 10)
    c1[1] ~ InverseGamma(2, 0.1)
    d[1] ~ InverseGamma(2, 0.0001)
    σ[1] ~ InverseGamma(2, 0.01)
    # Variance matrix
    Σ = σ[1].*I(N)
    # Find mean using the siler mortality function
    m_means = exp.(-b0[1] .- b1[1].* ages) .+ exp.(-c0[1] .+ c1[1] .* ages) .+ d[1]
    # Draw from truncated normal dist
    m_data[1] ~ MvNormal(m_means, Σ)
    # Time series innovations
    ϵ_b0 ~ filldist(Normal(0, σ_b0), T-1)
    ϵ_b1 ~ filldist(Normal(0, σ_b1), T-1)
    ϵ_c0 ~ filldist(Normal(0, σ_c0), T-1)
    ϵ_c1 ~ filldist(Normal(0, σ_c1), T-1)
    ϵ_d ~ filldist(Normal(0, σ_d), T-1)
    ϵ_σ ~ filldist(Normal(0, σ_σ), T-1)
    for tt in 2:T
        # Update parameters
        b0[tt] =  max(b0[tt-1] + ϵ_b0[tt-1], 0.0)
        b1[tt] = max(b1[tt-1] + ϵ_b1[tt-1], 0.0)
        c0[tt] = max(c0[tt-1] + ϵ_c0[tt-1], 0.0)
        c1[tt] = max(c1[tt-1] + ϵ_c1[tt-1], 0.0)
        d[tt] = max(d[tt-1] + ϵ_d[tt-1], 0.0)
        σ[tt] = max(σ[tt-1] + ϵ_σ[tt-1], 0.0)
        # Variance matrix
        Σ = σ[tt].*I(N)
        # Find mean using the siler mortality function
        m_means = exp.(-b0[tt] .- b1[tt].* ages) .+ exp.(-c0[tt] .+ c1[tt] .* ages) .+ d[tt]
        # Draw from truncated normal dist
        m_data[tt] ~ MvNormal(m_means, Σ)
    end
end


periods = 1:87
periods = Int.(1:5:86)
#periods = [1,10,20,30,40,50,60,70,80]
years_selected = years[periods]
iterations = 25
# MAP estimate to initialise MCMC
map_dyn = optimize(siler_dyn(m_data[periods], ages), MAP(),
    Optim.Options(iterations=1_000))
display(coef(map_dyn)[1:12])
display(map_dyn[Symbol("b0[1]")])
# Sample
chain_dyn = sample(siler_dyn(m_data[periods], ages), NUTS(0.65), iterations,
    init_params = map_dyn.values.array)
display(chain_dyn)

#chain_dyn1 = sample(siler_dyn(m_data[periods], ages), NUTS(0.65), iterations, chn=4)
#chain_dyn2 = sample(siler_dyn(m_data[periods], ages), NUTS(0.65), iterations, chn=4)
#chain_dyn3 = sample(siler_dyn(m_data[periods], ages), NUTS(0.65), iterations, chn=4)

"""
Calculate the implied latent parameter values from shocks
"""
df_dyn = DataFrame(chain_dyn)
for tt in 1:(length(periods)-1)
    # b0
    b0_temp = df_dyn[:,Symbol("b0["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_b0["*string(tt)*"]")]
    df_dyn[:, Symbol("b0["*string(tt+1)*"]")] = b0_temp
    # b1
    b1_temp = max.(df_dyn[:,Symbol("b1["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_b1["*string(tt)*"]")], [0.0])
    df_dyn[:, Symbol("b1["*string(tt+1)*"]")] = b1_temp
    # c0
    c0_temp = df_dyn[:,Symbol("c0["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_c0["*string(tt)*"]")]
    df_dyn[:, Symbol("c0["*string(tt+1)*"]")] = c0_temp
    # c1
    c1_temp = max.(df_dyn[:,Symbol("c1["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_c1["*string(tt)*"]")], [0.0])
    df_dyn[:, Symbol("c1["*string(tt+1)*"]")] = c1_temp
    # d
    d_temp = max.(df_dyn[:,Symbol("d["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_d["*string(tt)*"]")], [0.0])
    df_dyn[:, Symbol("d["*string(tt+1)*"]")] = d_temp
    # σ
    σ_temp = max.(df_dyn[:,Symbol("σ["*string(tt)*"]")] .+ df_dyn[:,Symbol("ϵ_σ["*string(tt)*"]")], [0.0])
    df_dyn[:, Symbol("σ["*string(tt+1)*"]")] = σ_temp
end

b0s = vcat([:iteration, :chain], Symbol.("b0[".*string.(1:length(periods)).*"]"))
b1s = vcat([:iteration, :chain], Symbol.("b1[".*string.(1:length(periods)).*"]"))
c0s = vcat([:iteration, :chain], Symbol.("c0[".*string.(1:length(periods)).*"]"))
c1s = vcat([:iteration, :chain], Symbol.("c1[".*string.(1:length(periods)).*"]"))
ds = vcat([:iteration, :chain], Symbol.("d[".*string.(1:length(periods)).*"]"))
σs = vcat([:iteration, :chain], Symbol.("σ[".*string.(1:length(periods)).*"]"))

# Extract the posterior distributions and some summary statistics
b0_ests = summarise_stats(df_dyn[:,b0s], years_selected)
b1_ests = summarise_stats(df_dyn[:,b1s], years_selected)
c0_ests = summarise_stats(df_dyn[:,c0s], years_selected)
c1_ests = summarise_stats(df_dyn[:,c1s], years_selected)
d_ests = summarise_stats(df_dyn[:,ds], years_selected)
σ_ests = summarise_stats(df_dyn[:,σs], years_selected)



plot(layout = (2,3))
plot!(b0_ests.variable, -b0_ests.median./b1_ests.median, title = L"B_{t}", label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, -b0_ests.pc975./b1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, -b0_ests.pc025./b1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, -b0_ests.pc85./b1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, -b0_ests.pc15./b1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)

plot!(b1_ests.variable, b1_ests.median, title = L"b_{t}", label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)

plot!(c0_ests.variable, -c0_ests.median./c1_ests.median, title = L"C_{t}", label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, -c0_ests.pc975./c1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, -c0_ests.pc025./c1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, -c0_ests.pc85./c1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, -c0_ests.pc15./c1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 4)

plot!(c1_ests.variable, c1_ests.median, title = L"c_{t}", label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 5)

plot!(d_ests.variable, d_ests.median, title = L"d_{t}", label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 3)

plot!(σ_ests.variable, σ_ests.median, title = L"\sigma_{t}", label = false, color = 1, subplot = 6)#, ylim = (0.0, 0.1))
plot!(σ_ests.variable, σ_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 6)

savefig("figures/Siler_dynamic/siler_dyn.pdf")




# Plot model fit
plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
scatter!(m_data[1], markershape = :cross, markeralpha = 0.5, label = "Data (1933)")
plot!(siler(median(chain_1933[:b0]), median(chain_1933[:b1]), median(chain_1933[:c0]),
    median(chain_1933[:c1]), median(chain_1933[:d]), ages), label = "Siler MCMC fit (1993)")
scatter!(m_data[87], markershape = :xcross, markeralpha = 0.5, label = "Data (2019)")
plot!(siler(median(chain_2019[:b0]), median(chain_2019[:b1]), median(chain_2019[:c0]),
    median(chain_2019[:c1]), median(chain_2019[:d]), ages), label = "Siler MCMC fit (2019)")
savefig("figures/Siler_static/siler_fit.pdf")






"""
Multiple static Siler models
"""
# Define Turing model
@model function siler_static_T(m_data, ages)
    T = length(m_data)

    # Our prior beliefs
    b0 ~ filldist(InverseGamma(2, 10), T)
    b1 ~ filldist(InverseGamma(2, 0.1), T)
    c0 ~ filldist(InverseGamma(2, 10), T)
    c1 ~ filldist(InverseGamma(2, 0.1), T)
    d ~ filldist(InverseGamma(2, 0.001), T)
    σ ~ filldist(InverseGamma(2, 0.001), T)

    # The number of observations.
    for tt in 1:T
        Σ = σ[tt].*I(N)
        # Find mean using the siler mortality function
        m_means = exp.(-b0[tt] .- b1[tt].* ages) .+ exp.(-c0[tt] .+ c1[tt] .* ages) .+ d[tt]
        # Draw from truncated normal dist
        m_data[tt] ~ MvNormal(m_means, Σ)
        # Find mean using the siler mortality function
    end
end

# Estimate the model
periods = Int.(1:5:86)
years_selected = years[periods]
iterations = 5000
chain_staticT = sample(siler_static_T(m_data[periods]), NUTS(0.65), iterations,chn=4)
display(chain_staticT)


prior_staticT = sample(siler_static_T(m_data[periods]), Prior(), 100)

# Extract the posterior distributions and some summary statistics
df_staticT = DataFrame(chain_staticT)
b0_ests = summarise_stats(df_staticT[:,b0s], years_selected)
b1_ests = summarise_stats(df_staticT[:,b1s], years_selected)
c0_ests = summarise_stats(df_staticT[:,c0s], years_selected)
c1_ests = summarise_stats(df_staticT[:,c1s], years_selected)
d_ests = summarise_stats(df_staticT[:,ds], years_selected)
σ_ests = summarise_stats(df_staticT[:,σs], years_selected)


plot(layout = (2,3))
plot!(b0_ests.variable, b0_ests.median, title = L"b_{0}", label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
plot!(b0_ests.variable, b0_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)

plot!(b1_ests.variable, b1_ests.median, title = L"b_{1}", label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
plot!(b1_ests.variable, b1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)

plot!(c0_ests.variable, c0_ests.median, title = L"c_{0}", label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 4)
plot!(c0_ests.variable, c0_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 4)

plot!(c1_ests.variable, c1_ests.median, title = L"c_{1}", label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 5)
plot!(c1_ests.variable, c1_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 5)

plot!(d_ests.variable, d_ests.median, title = L"d", label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 3)
plot!(d_ests.variable, d_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 3)

plot!(σ_ests.variable, σ_ests.median, title = L"\sigma", label = false, color = 1, subplot = 6)#, ylim = (0.0, 0.1))
plot!(σ_ests.variable, σ_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 6)
plot!(σ_ests.variable, σ_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 6)

savefig("figures/Siler_dynamic/siler_staticT.pdf")

















##################################################################################










df = CSV.read("data/data_master_1.csv", DataFrame, ntasks = false)
df[:,:sp_500]
s = df.sp_500
s_diff = diff(s)

@model function ARIMA011(x)
    T = length(x)
    # Set up error vector.
    ϵ = Vector(undef, T)
    x_hat = Vector(undef, T)
    θ ~ Uniform(-5, 5)
    # Treat the first x_hat as a parameter to estimate.
    x_hat[1] ~ Normal(0, 1)
    ϵ[1] = x[1] - x_hat[1]
    for t in 2:T
        # Predicted value for x.
        x_hat[t] = x[t-1] - θ * ϵ[t-1]
        # Calculate observed error.
        ϵ[t] = x[t] - x_hat[t]
        # Observe likelihood.
        x[t] ~ Normal(x_hat[t], 1)
    end
end


chain_ARIMA011 = sample(ARIMA011(s), NUTS(0.6), 1000)





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
