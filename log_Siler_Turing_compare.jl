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

"""
Import and clean data
"""
## Import the mortality data
all_df = CSV.read("data/clean/all_lifetab.csv", DataFrame, ntasks = 1)

# Best practice data
mort_df = all_df[(all_df.year .>= 1900), :]
sort!(mort_df, [:code, :year, :age])
bp_df = mort_df[(mort_df.best_practice .== 1), :]
plot(bp_df.age, bp_df.mx, group = bp_df.year, legend = :top)

# Restrict to G7 post 1900 for now
G7_countries = ["Canada", "France", "West Germany", "Italy", "Japan", "United Kingdom",
    "United States of America"]
G7_df = mort_df[in.(mort_df.name, [G7_countries]),:]


"""
Preliminary checks and illustrative estimation of first and last period
"""
## Data prep for single coujntry
country_df = G7_df[(G7_df.code .== "USA"), :]
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


## Estimate independent models for first and last periods
periods = Int.([1,T])
years_selected = Int.(round.(country_years[periods]))
iterations = 2000
# Find MAP estimate as starting point
@time map_indep = optimize(log_siler_indep(country_lm_data[periods], country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
# MCMC sampling
@time chain_indep = sample(log_siler_indep(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    iterations, 4, init_params = map_indep.values.array)
display(chain_indep)
# Plot sampled posteriors across threads
plot(chain_indep[["B[1]", "b[1]", "C[1]", "c[1]", "d[1]", "σ[1]"]])
plot(chain_indep[["B[2]", "b[2]", "C[2]", "c[2]", "d[2]", "σ[2]"]])
# Plot parameters as time series
parests_indep = extract_variables(chain_indep, years_selected)
plot_siler_params(parests_indep)
# Visualise model fit
plot_fit_year(parests_indep, country_m_data[1], years_selected[1], log_vals = false)
plot_fit_year!(parests_indep, country_m_data[T], years_selected[end], log_vals = false)



"""
Estimate Siler model for each country in G7
"""
## Set some estimation options
niters = 2000
nthreads = 4

parests_all = DataFrame(code = Symbol[], name = String[], parameter = Symbol[], year = Int64[],
    mean = Float64[], min = Float64[], median = Float64[], max = Float64[],
    nmissing  = Int64[], eltype = DataType[], std = Float64[],
    pc975 = Float64[], pc025 = Float64[], pc85 = Float64[], pc15 = Float64[],
    pc75 = Float64[], pc25 = Float64[], mcse = Float64[], ess = Float64[], rhat = Float64[],
    ess_per_sec = Float64[])

parests_dict = Dict{Symbol, DataFrame}()

for code in unique(G7_df.code)
    print("Working on model for "*code)
    # Extract and convert relevant data into correct form
    country_df = G7_df[(G7_df.code .== code), :]
    country_m_data = chunk(country_df.mx, 110)
    country_lm_data = [log.(m_dist) for m_dist in country_m_data]
    country_ages = Int64.(0:maximum(country_df.age))
    country_years = unique(country_df.year)
    T = length(country_lm_data)
    name = country_df.name[1]

    @assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
    @assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"

    # MAP estimate to initialise MCMC
    @time map_dyn = optimize(log_siler_dyn(country_lm_data, country_ages), MAP(), LBFGS(),
        Optim.Options(iterations=50_000, allow_f_increases=true))
    print("Estimated MAP for "*code)
    # Estimate by MCMC
    @time chain_dyn = sample(log_siler_dyn(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
        niters, nthreads, init_params = map_dyn.values.array)
    print("Sampled posterior for "*code)
    # Plot some example posterior distributions
    display(chain_dyn)
    plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]", "lc[1]", "ld[1]", "lσ[1]"]])
    plot!(margin=8Plots.mm)
    savefig("results/"*code*"/"*code*"_first_period_posteriors.pdf")
    plot(chain_dyn[["σ_pars[1]", "σ_pars[2]", "σ_pars[3]", "σ_pars[4]", "σ_pars[5]", "σ_pars[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results/"*code*"/"*code*"_rw_variance_posteriors.pdf")
    # Extract and plot results
    parests_dyn = extract_variables(chain_dyn, country_years, log_pars = true,
        σ_pars = true)
    p1 = plot_siler_params(parests_dyn)
    p_title = plot(title = "Dynamic Siler parameters "*string(code), grid = false, showaxis = false,
        bottom_margin = -10Plots.px, yticks = false, xticks = false)
    display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
    savefig("results/"*code*"/"*code*"_param_estimates.pdf")

    plot_fit_year(parests_dyn, country_m_data[1], country_years[1], log_vals = false, col = 1)
    plot_fit_year!(parests_dyn, country_m_data[T], country_years[T], log_vals = false, col = 2)
    plot!(title = "Siler model fit "*string(code))
    savefig("results/"*code*"/"*code*"_model_fit.pdf")

    # Store and export results
    CSV.write("results/"*code*"/"*code*"_est_results.csv", parests_dyn)
    insertcols!(parests_dyn, 1, :code => repeat([Symbol(code)], nrow(parests_dyn)) )
    insertcols!(parests_dyn, 2, :name => repeat([name], nrow(parests_dyn)) )

    parests_all = vcat(parests_all, parests_dyn)
    sort!(parests_all, [:code])
    CSV.write("results/country_siler_est_results.csv", parests_all)

    parests_dict[Symbol(code)] = parests_dyn

end



"""
Estimate Siler model for best practice as a comparison
"""
# Sort bp_df by year and age
code = "BestPractice"
name = "Best Practice"
sort!(bp_df, [:year, :age])
# Extract and convert relevant data into correct form
country_m_data = chunk(bp_df.mx, 110)
country_lm_data = [log.(m_dist) for m_dist in country_m_data]
country_ages = Int64.(0:maximum(bp_df.age))
country_years = unique(bp_df.year)
T = length(country_lm_data)

@assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
@assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"

# MAP estimate to initialise MCMC
@time map_dyn = optimize(log_siler_dyn(country_lm_data, country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
print("Estimated MAP for "*code)
# Estimate by MCMC
@time chain_dyn = sample(log_siler_dyn(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
    niters, nthreads, init_params = map_dyn.values.array)
print("Sampled posterior for "*code)
# Plot some example posterior distributions
display(chain_dyn)
plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]", "lc[1]", "ld[1]", "lσ[1]"]])
plot!(margin=8Plots.mm)
savefig("results/"*code*"/"*code*"_first_period_posteriors.pdf")
plot(chain_dyn[["σ_pars[1]", "σ_pars[2]", "σ_pars[3]", "σ_pars[4]", "σ_pars[5]", "σ_pars[6]"]])
plot!(margin=8Plots.mm)
savefig("results/"*code*"/"*code*"_rw_variance_posteriors.pdf")
# Extract and plot results
parests_dyn = extract_variables(chain_dyn, country_years, log_pars = true,
    σ_pars = true)
p1 = plot_siler_params(parests_dyn)
p_title = plot(title = "Dynamic Siler parameters "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("results/"*code*"/"*code*"_param_estimates.pdf")

plot_fit_year(parests_dyn, country_m_data[1], country_years[1], log_vals = false, col = 1)
plot_fit_year!(parests_dyn, country_m_data[T], country_years[T], log_vals = false, col = 2)
plot!(title = "Siler model fit "*string(code))
savefig("results/"*code*"/"*code*"_model_fit.pdf")

# Store and export results
CSV.write("results/"*code*"/"*code*"_est_results.csv", parests_dyn)
insertcols!(parests_dyn, 1, :code => repeat([Symbol(code)], nrow(parests_dyn)) )

parests_all = vcat(parests_all, parests_dyn)
sort!(parests_all, [:code])
CSV.write("results/country_siler_est_results.csv", parests_all)

parests_dict[Symbol(code)] = parests_dyn


"""
Merge the names back in (should have done this at the beginning)
"""
insertcols!(parests_all, 2, :name => repeat([" "], nrow(parests_all)) )

parests_all.name[parests_all.code .== :BestPractice] .= "Best Practice"
parests_all.name[parests_all.code .== :CAN] .= "Canada"
parests_all.name[parests_all.code .== :DEUTW] .= "West Germany"
parests_all.name[parests_all.code .== :FRA] .= "France"
parests_all.name[parests_all.code .== :GBR] .= "United Kingdom"
parests_all.name[parests_all.code .== :ITA] .= "Italy"
parests_all.name[parests_all.code .== :JPN] .= "Japan"
parests_all.name[parests_all.code .== :USA] .= "United States of America"
