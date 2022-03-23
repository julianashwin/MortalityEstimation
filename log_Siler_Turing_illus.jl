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
# Choose year
yy = 1
## Find some starting points
# MAP estimate for static model on raw mortality
map_static = optimize(siler_static((country_m_data[yy]), country_ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
map_vals =  map_static.values.array
map_param = SilerParam(b = map_vals[2], B = map_vals[1], c = map_vals[4], C = map_vals[3],
    d = map_vals[5], σ = map_vals[6])
# Map estimate for static model on log mortality
map_log_static = optimize(log_siler_static((country_lm_data[yy]), country_ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
map_log_vals =  map_log_static.values.array
map_log_param = SilerParam(b = map_log_vals[2], B = map_log_vals[1], c = map_log_vals[4], C = map_log_vals[3],
    d = map_log_vals[5], σ = map_log_vals[6])
# Alternatively, we can simulate from the prior and start there
prior_static = sample(siler_static(country_m_data[yy], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_static)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_vals = df_prior[1,3:(end-1)]
prior_param = SilerParam(b = prior_vals.b, B = prior_vals.B, c = prior_vals.c, C = prior_vals.C,
    d = prior_vals.d, σ = prior_vals.σ)


## Estimate models
# Number of MCMC iterations and chains
niters = 800
nchains = 1
# Raw mortality model
chain_static = sample(siler_static(country_m_data[1], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = map_log_vals)
df_chain = DataFrame(chain_static)
insert!.(eachcol(df_chain), 1, vcat([0,0],median.(eachcol(df_chain[:,3:end]))))
chain_vals = df_chain[1,:]
chain_param = SilerParam(b = chain_vals.b, B = chain_vals.B, c = chain_vals.c, C = chain_vals.C,
    d = chain_vals.d, σ = chain_vals.σ)
# Raw mortality model
chain_log_static = sample(log_siler_static(country_lm_data[1], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = map_log_vals)
df_chain = DataFrame(chain_log_static)
insert!.(eachcol(df_chain), 1, vcat([0,0],median.(eachcol(df_chain[:,3:end]))))
chain_log_vals = df_chain[1,:]
chain_log_param = SilerParam(b = chain_log_vals.b, B = chain_log_vals.B, c = chain_log_vals.c, C = chain_log_vals.C,
    d = chain_log_vals.d, σ = chain_log_vals.σ)

## Plot fit
scatter(country_m_data[yy], markershape = :cross, markeralpha = 0.5,
    label = "Data ("*string(country_years[1])*")", legend = :topleft)
plot!(siler.([prior_param], 0:110, spec = :Colchero), label = "Static model prior")
plot!(siler.([map_param], 0:110, spec = :Colchero), label = "Static model MAP", linestyle = :dash)
plot!(siler.([chain_param], 0:110, spec = :Colchero), label = "Static model post. median", linestyle = :dot)
plot!(siler.([map_log_param], 0:110, spec = :Colchero), label = "Static log model MAP", linestyle = :dash)
plot!(siler.([chain_log_param], 0:110, spec = :Colchero), label = "Static log model post. median", linestyle = :dot)


"""
Multiple independent Siler models
"""
# Adjust if you don't want every period
periods = Int.(1:7:T)
years_selected = Int.(round.(country_years[periods]))

## Find some starting points
# MAP estimate for multiple independent models on log mortality
map_indep = optimize(log_siler_indep(country_lm_data[periods], country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=60_000, allow_f_increases=true))
map_indep_vals =  map_indep.values.array
# Alternatively, we can simulate from the prior and start there
prior_indep = sample(log_siler_indep(country_lm_data[periods], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_indep)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_indep_vals = df_prior[1,3:(end-1)]

## Estimate the model
niters = 800
nchains = 1
# MCMC sampling
chain_indep = sample(log_siler_indep(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_indep_vals)
display(chain_indep)
parests_indep = extract_variables(chain_indep, years_selected, spec = :Colchero)
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
