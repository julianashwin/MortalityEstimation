"""
Use MCMC from Turing.jl to estimate log Siler model on log mortality data
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
using TruncatedDistributions, PDMats
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings

include("src/MortalityEstimation.jl")


## Import the mortality data
all_df = CSV.read("data/clean/all_lifetab.csv", DataFrame, ntasks = 1)

# Extract the best practice series
bp_df = all_df[(all_df.best_practice .== 1) .& (all_df.year .>= 1900), :]
sort!(bp_df, [:year, :age])
bp_alt_df = all_df[(all_df.best_practice_alt .== 1) .& (all_df.year .>= 1900), :]
sort!(bp_alt_df, [:year, :age])
# Also a special case for whole sample period in SWE
swe_df = all_df[(all_df.code .== "SWE"),:]
sort!(swe_df, [:year, :age])
# And a general post 1900 df
mort_df = all_df[(all_df.year .>= 1900), :]
sort!(mort_df, [:code, :year, :age])
showtable(mort_df)


"""
Plot priors
"""
## Set priors
# Mean of IG is b/(a-1)
plot(layout = (2,3), yticks = false, size = (1000,500))
plot!(LogNormal(log(10), 2.0), title = L"B_{1} \sim LogNormal(ln(10),2)",
    label = false, subplot = 1, xlims = (0,100))
plot!(LogNormal(log(2), 1.0), xlim=(0,10), title = L"b_{1} \sim LogNormal(ln(2),2)",
    label = false, subplot = 2)
plot!(LogNormal(log(10), 2.0), title = L"C_{1} \sim LogNormal(ln(120),2)",
    label = false, subplot = 4, xlims = (0,200))
plot!(LogNormal(log(0.1), 1.0), xlim=(0,1), title = L"c_{1} \sim LogNormal(ln(0.1),2)",
    label = false, subplot = 5)
plot!(LogNormal(log(0.025), 1.0), title = L"d_{1} \sim LogNormal(ln(0.025),2)",
    label = false, subplot = 3)
plot!(LogNormal(log(0.001), 1.0), title = L"\sigma_{1} \sim LogNormal(ln(0.001),2)",
    label = false, subplot = 6)
savefig("figures/general/log_siler_priors.pdf")


"""
Set up a single country/case to run through examples
"""
## Data prep for single coujntry
code = "SWE"
country_df = mort_df[(mort_df.code .== code), :]
country_df = bp_df
# Check data looks sensible
plot(country_df.age, country_df.mx, group = country_df.year, legend = :top)
# Convert this into a matrix of mortality rates over time, age and year vectors
country_m_data = chunk(country_df.mx, 110)
country_lm_data = [log.(m_dist) for m_dist in country_m_data]
country_ages = Int64.(0:maximum(country_df.age))
country_years = unique(country_df.year)
T = length(country_lm_data)

@assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
@assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"



"""
Static Siler model
"""
## Declare our Turing model for Siler
@model function log_siler_static(lm_dist, ages)
    # The number of observations.
    N = length(lm_dist)
    # Our prior beliefs
    B ~ LogNormal(log(10), 2.0)
    b ~ LogNormal(log(2), 1.0)
    C ~ LogNormal(log(120), 2.0)
    c ~ LogNormal(log(0.1), 1.0)
    d ~ LogNormal(log(0.025), 1.0)
    σ ~ LogNormal(log(0.001), 1.0)
    # Define the logged parameters, which should be normally distributed
    lB = log(B)
    lb = log(b)
    lC = log(C)
    lc = log(c)
    ld = log(d)
    lσ = log(σ)
    # Find mean using the siler mortality function
    μs = exp.(- exp(lb).* (ages .+ exp(lB))) .+ exp.(exp(lc) .* (ages .- exp(lC))) .+ exp(ld)
    m_vars = exp(σ).*ones(N)
    #m_vars[m_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ = Diagonal(m_vars)
    # Draw from normal dist
    lm_dist ~ MvNormal(log.(μs), Σ)
    #m_dist = exp.(lm_dist)
end

## Estimate models
# Number of MCMC iterations
iterations = 2500
# Sample for first period
map_static_1 = optimize(log_siler_static(log.(m_data[1]), ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
@time chain_1 = sample(log_siler_static(log.(m_data[1]), ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_static_1.values.array)
display(chain_1)
plot(chain_1)
# Sample for last period
map_static_T = optimize(log_siler_static(log.(m_data[T]), ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
@time chain_T = sample(log_siler_static(log.(m_data[T]), ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_static_T.values.array)
display(chain_T)
plot(chain_T)

## Plot model fit
plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
scatter!(m_data[1], markershape = :cross, markeralpha = 0.5, label = "Data ("*string(years[1])*")")
plot!(siler(median(chain_1[:B]), median(chain_1[:b]), median(chain_1[:C]),
    median(chain_1[:c]), median(chain_1[:d]), ages), label = "Siler MCMC fit ("*string(years[1])*")")
scatter!(m_data[T], markershape = :xcross, markeralpha = 0.5, label = "Data ("*string(years[T])*")")
plot!(siler(median(chain_T[:B]), median(chain_T[:b]), median(chain_T[:C]),
    median(chain_T[:c]), median(chain_T[:d]), ages), label = "Siler MCMC fit ("*string(years[T])*")")
savefig("figures/Siler_static/log_siler_fit.pdf")

plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Log Mortality")
scatter!(log.(m_data[1]), markershape = :cross, markeralpha = 0.5, label = "Data ("*string(years[1])*")")
plot!(log.(siler(median(chain_1[:B]), median(chain_1[:b]), median(chain_1[:C]),
    median(chain_1[:c]), median(chain_1[:d]), ages)),
    label = "Siler MCMC fit ("*string(years[1])*")")
scatter!(log.(m_data[T]), markershape = :xcross, markeralpha = 0.5, label = "Data ("*string(years[T])*")")
plot!(log.(siler(median(chain_T[:B]), median(chain_T[:b]), median(chain_T[:C]),
    median(chain_T[:c]), median(chain_T[:d]), ages)),
    label = "Siler MCMC fit ("*string(years[T])*")")
savefig("figures/Siler_static/log_siler_logfit.pdf")



## Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
plot(layout = (2,3), size = (800, 400))
density!(vcat(chain_1[:B]...), title = L"B", label = string(years[1]), subplot = 1)
density!(vcat(chain_T[:B]...), label = string(years[T]), subplot = 1, legend = :topright)
density!(vcat(chain_1[:b]...), title = L"b", subplot = 2, legend = false)
density!(vcat(chain_T[:b]...), subplot = 2, legend = false)
density!(vcat(chain_1[:C]...), title = L"C", subplot = 4, legend = false)
density!(vcat(chain_T[:C]...), subplot = 4, legend = false)
density!(vcat(chain_1[:c]...), title = L"c", subplot = 5, legend = false)
density!(vcat(chain_T[:c]...), subplot = 5, legend = false)
density!(vcat(chain_1[:d]...), title = L"d", subplot = 3, legend = false)
density!(vcat(chain_T[:d]...), subplot = 3, legend = false)
density!(vcat(chain_1[:σ]...), title = L"\sigma", subplot = 6, legend = false)
density!(vcat(chain_T[:σ]...), subplot = 6, legend = false, xrotation = 45.0, margin=3Plots.mm)
savefig("figures/Siler_static/log_siler_1vsT.pdf")


"""
Multiple independent Siler models
"""
## Define Turing model
@model function log_siler_indep(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Our prior beliefs
    B ~ filldist(LogNormal(log(10), 2.0), T)
    b ~ filldist(LogNormal(log(2), 1.0), T)
    C ~ filldist(LogNormal(log(120), 2.0), T)
    c ~ filldist(LogNormal(log(0.1), 1.0), T)
    d ~ filldist(LogNormal(log(0.025), 1.0), T)
    σ ~ filldist(LogNormal(log(0.1), 1.0), T)
    # Define the logged parameters, which should be normally distributed
    lB = log.(B)
    lb = log.(b)
    lC = log.(C)
    lc = log.(c)
    ld = log.(d)
    lσ = log.(σ)

    for tt in 1:T
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end

## Estimate the model
periods = Int.(1:T)
years_selected = Int.(round.(years[periods]))
iterations = 4000
# Find MAP estimate as starting point
lm_data = [log.(m_dist) for m_dist in m_data]
@time map_indep = optimize(log_siler_indep(lm_data[periods], ages), MAP(), LBFGS(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
# MCMC sampling
@time chain_indep = sample(log_siler_indep(lm_data[periods], ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_indep.values.array)
display(chain_indep)

parests_indep = extract_variables(chain_indep, periods, years_selected)

plot_siler_params(parests_indep)
savefig("figures/Siler_staticT/log_siler_staticT.pdf")


## Visualise model fit

function plot_fit_year(parests, m_dist, year; log_vals = false)

    # Extract the parameters for that year
    B = parests.mean[(parests.year .== year).*(parests.parameter .== :B)][1]
    b = parests.mean[(parests.year .== year).*(parests.parameter .== :b)][1]
    C = parests.mean[(parests.year .== year).*(parests.parameter .== :C)][1]
    c = parests.mean[(parests.year .== year).*(parests.parameter .== :c)][1]
    d = parests.mean[(parests.year .== year).*(parests.parameter .== :d)][1]
    σ = parests.mean[(parests.year .== year).*(parests.parameter .== :σ)][1]

    plt = plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
    if log_vals
        scatter!(log.(m_dist), markershape = :cross, markeralpha = 0.5,
            label = "Data ("*string(year)*")", color = :black)
        plot!(log.(siler(B,b,C,c,d, ages)),
            label = "Siler MCMC fit ("*string(year)*")", color = :blue)
        plot!(log.(siler(B,b,C,c,d, ages)) .+ 2*σ,
            label = L"\pm\ 2\ s.e.", color = :blue, linestyle = :dash)
        plot!(log.(siler(B,b,C,c,d, ages)) .- 2*σ,
            label = false, color = :blue, linestyle = :dash)
    else
        scatter!((m_dist), markershape = :cross, markeralpha = 0.5, label = "Data ("*string(year)*")")
        plot!((siler(B,b,C,c,d, ages)), label = "Siler MCMC fit ("*string(year)*")")
    end

    return plt
end

plot_fit_year(parests, m_data[1], Int.(round.(years_selected))[1], log_vals = true)

plot_fit_year(parests, m_data[T-3], Int.(round.(years_selected))[T-3], log_vals = true)


"""
Dynamic Siler model
"""
# Declare our Turing model for dynamic Siler equation
@model function log_siler_dyn(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Parameters
    B = Vector(undef, T)
    b = Vector(undef, T)
    C = Vector(undef, T)
    c = Vector(undef, T)
    d = Vector(undef, T)
    σ = Vector(undef, T)
    # Logged parameters for random walk model
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.1),6)
    # Priors on drift terms
    #μ_pars ~ filldist(Normal(0, 0.1), 6)

    # First period from priors
    lB[1] ~ Normal(log(10), 2.0)
    lb[1] ~ Normal(log(2), 1.0)
    lC[1] ~ Normal(log(120), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.025), 1.0)
    lσ[1] ~ Normal(log(0.1), 1.0)
    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*(ages .+ exp(lB[1]))) .+
        exp.(exp(lc[1]).*(ages .- exp(lC[1]))) .+ exp(ld[1])
    lm_vars = exp(lσ[1]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ)
    # Loop through random walk process
    for tt in 2:T
        # Calculate updated variances
        #Σ_pars = Diagonal(σ_pars)
        #lmean_pars = [lB[tt-1], lb[tt-1], lC[tt-1], lc[tt-1], ld[tt-1], lσ[tt-1]]
        var_B = max(σ_pars[1], 1e-8)
        var_b = max(σ_pars[2], 1e-8)
        var_C = max(σ_pars[3], 1e-8)
        var_c = max(σ_pars[4], 1e-8)
        var_d = max(σ_pars[5], 1e-8)
        var_σ = max(σ_pars[6], 1e-8)
        # Update parameters
        lB[tt] ~ Normal(lB[tt-1], var_B)
        lb[tt] ~ Normal(lb[tt-1], var_b)
        lC[tt] ~ Normal(lC[tt-1], var_C)
        lc[tt] ~ Normal(lc[tt-1], var_c)
        ld[tt] ~ Normal(ld[tt-1], var_d)
        lσ[tt] ~ Normal(lσ[tt-1], var_σ)
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end


periods = Int.(1:T)
#periods = Int.((T-30):10:T)
#periods = [1,10,20,30,40,50,60,70,80]
years_selected = Int.(round.(years[periods]))
iterations = 2000
lm_data = [log.(m_dist) for m_dist in m_data]

# MAP estimate to initialise MCMC
@time map_dyn = optimize(log_siler_dyn(lm_data[periods], ages), MAP(), LBFGS(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
display(exp.(coef(map_dyn)[80:110]))
# Estimate by MCMC
@time chain_dyn = sample(log_siler_dyn(lm_data[periods], ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_dyn.values.array)
display(chain_dyn)
plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]", "lc[1]", "ld[1]", "lσ[1]"]])
plot(chain_dyn[["σ_pars[1]", "σ_pars[2]", "σ_pars[3]", "σ_pars[4]", "σ_pars[5]", "σ_pars[6]"]])
plot!(margin=8Plots.mm)
savefig("figures/Siler_dynamic/rw_variance_posteriors.pdf")
#plot(chain_dyn[[:μ_B, :μ_b, :μ_C, :μ_c, :μ_d, :μ_σ]])
parests_dyn = extract_variables(chain_dyn, periods, years_selected, log_pars = true)

# Plot each parameter over time
plot_siler_params(parests_dyn)
savefig("figures/Siler_dynamic/siler_dyn.pdf")


plot_fit_year(parests_dyn, m_data[1], years_selected[1], log_vals = false)

plot_fit_year(parests_dyn, m_data[T-3], years_selected[T-3], log_vals = false)

savefig("figures/Siler_dynamic/siler_fit.pdf")
