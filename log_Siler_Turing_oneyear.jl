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
all_df = CSV.read("data/clean/bp.csv", DataFrame, ntasks = 1)

# Best practice data
mort_df = all_df[(all_df.year .>= 1900), :]
sort!(mort_df, [:code, :year, :age])
bp_df = mort_df #[(mort_df.best_practice .== 1), :]
plot(bp_df.age, bp_df.mx, group = bp_df.year, legend = :top)
sort!(bp_df, [:year, :age])



"""
Preliminary checks and illustrative estimation of first and last period
"""
## Data prep for single coujntry
code = "Best Practice"
country_df = bp_df
# Check data looks sensible
plot(country_df.age, country_df.mx, group = country_df.year, legend = :top)
# Suprisingly, we actually have some zeros here for small countries (e.g. ISL)
country_df.mx_f[country_df.mx_f .== 0.0] .=  minimum(country_df.mx_f[country_df.mx_f .> 0.0])
# Convert this into a matrix of mortality rates over time, age and year vectors
country_m_data = chunk(country_df.mx_f, 110)
country_lm_data = [log.(m_dist) for m_dist in country_m_data]
country_ages = Int64.(0:maximum(country_df.age))
country_years = unique(country_df.year)
T = length(country_lm_data)

@assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
@assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"


prior_est = sample(log_siler_dyn_i2drift(country_lm_data, country_ages), Prior(), 5000)
df_prior = DataFrame(prior_est)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_vals = df_prior[1,3:(end-1)]


# Estimate by MCMC
@time chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
    1000, 4, init_params = prior_vals)
display(chain_i2)
@save "figures/general/siler_i2_chain_oneyear.jld2" chain_i2

years_selected = country_years
# Display the parameters with Colchero specification
parests_i2_col = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
p1 = plot_siler_params(parests_i2_col)
p_title = plot(title = "I(2) Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/general/siler_i2_params_col.pdf")
CSV.write("figures/general/siler_i2_params_col.csv", parests_i2_col)
# With Scott specification
parests_i2_sco = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Scott)
p1 = plot_siler_params(parests_i2_sco)
p_title = plot(title = "I(2) Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/general/siler_i2_params_sco.pdf")
CSV.write("figures/general/siler_i2_params_sco.csv", parests_i2_sco)
# With Bergeron specification
parests_i2_ber = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Bergeron)
p1 = plot_siler_params(parests_i2_ber)
p_title = plot(title = "I(2) Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/general/siler_i2_params_ber.pdf")
CSV.write("figures/general/siler_i2_params_ber.csv", parests_i2_ber)
# With Standard specification
parests_i2_sta = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Standard)
p1 = plot_siler_params(parests_i2_sta)
p_title = plot(title = "I(2) Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/general/siler_i2_params_sta.pdf")
CSV.write("figures/general/siler_i2_params_sta.csv", parests_i2_sta)

## Plot time series parameters
p2 = plot_ts_params(parests_i2_col, model_vers = :i2drift)
p_title = plot(title = "I(2) Siler Colcher ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/general/siler_i2_ts_params_col.pdf")





"""
Estimate Siler model for each country in select
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

for code in unique(select_df.code)
    print("Working on model for "*code)
    # Extract and convert relevant data into correct form
    country_df = select_df[(select_df.code .== code), :]
    country_m_data = chunk(country_df.mx, 110)
    country_lm_data = [log.(m_dist) for m_dist in country_m_data]
    country_ages = Int64.(0:maximum(country_df.age))
    country_years = unique(country_df.year)
    T = length(country_lm_data)
    name = country_df.name[1]

    @assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
    @assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"

    # Sample from prior as a sense check
    chain_prior = sample(log_siler_dyn_firstdiff(country_lm_data, country_ages), Prior(), 5000)
    df_prior = DataFrame(chain_prior)
    insert!.(eachcol(df_prior), 1, vcat([0,0],mean.(eachcol(df_prior[:,3:end]))))
    prior_means = mean.(eachcol(df_prior[:,3:end]))

    # MAP estimate to initialise MCMC
    @time map_dyn = optimize(log_siler_dyn_firstdiff(country_lm_data, country_ages), MAP(), LBFGS(),
        Optim.Options(iterations=60_000, allow_f_increases=true))
    print("Estimated MAP for "*code)
    # Estimate by MCMC
    @time chain_dyn = sample(log_siler_dyn_firstdiff(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
        niters, nthreads, init_params = map_dyn.values.array)
    #@time chain_dyn = sample(log_siler_dyn_ext(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
    #    niters, nthreads, init_params =  Array(df_prior[1,3:(end-1)]))
    print("Sampled posterior for "*code)
    # Plot some example posterior distributions
    display(chain_dyn)
    plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]",
        "lc[1]", "ld[1]", "lσ[1]"]])
    plot!(margin=8Plots.mm)
    savefig("results_firstdiff/"*code*"/"*code*"_first_period_posteriors.pdf")
    plot(chain_dyn[["σ_par[1]", "σ_par[2]", "σ_par[3]", "σ_par[4]", "σ_par[5]", "σ_par[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results_firstdiff/"*code*"/"*code*"_rw_variance_posteriors.pdf")
    plot(chain_dyn[["α_pars[1]", "α_pars[2]", "α_pars[3]", "α_pars[4]", "α_pars[5]", "α_pars[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results_firstdiff/"*code*"/"*code*"_rw_drift_posteriors.pdf")
    plot(chain_dyn[["β_pars[1]", "β_pars[2]", "β_pars[3]", "β_pars[4]", "β_pars[5]", "β_pars[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results_firstdiff/"*code*"/"*code*"_rw_firstdiff_posteriors.pdf")

    # Extract and plot results
    parests_dyn = extract_variables(chain_dyn, country_years, log_pars = true,
        σ_pars = true, ext = true)
    p1 = plot_siler_params(parests_dyn)
    p_title = plot(title = "Dynamic (first diff) Siler parameters "*string(code), grid = false, showaxis = false,
        bottom_margin = -10Plots.px, yticks = false, xticks = false)
    display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
    savefig("results_firstdiff/"*code*"/"*code*"_param_estimates.pdf")
    # Plot the time series parameters
    p2 = plot_ts_params(parests_dyn)
    p_title = plot(title = "Siler time series (first diff) parameters "*string(code), grid = false, showaxis = false,
        bottom_margin = -10Plots.px, yticks = false, xticks = false)
    display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
    savefig("results_firstdiff/"*code*"/"*code*"_ts_params.pdf")

    plot_fit_year(parests_dyn, country_m_data[1], country_years[1], log_vals = false, col = 1)
    plot_fit_year!(parests_dyn, country_m_data[T], country_years[T], log_vals = false, col = 2)
    plot!(title = "Siler model (first diff) fit "*string(code))
    savefig("results_firstdiff/"*code*"/"*code*"_model_fit.pdf")

    # Store and export results
    CSV.write("results_firstdiff/"*code*"/"*code*"_est_results.csv", parests_dyn)
    insertcols!(parests_dyn, 1, :code => repeat([Symbol(code)], nrow(parests_dyn)) )
    insertcols!(parests_dyn, 2, :name => repeat([name], nrow(parests_dyn)) )

    parests_all = vcat(parests_all, parests_dyn)
    sort!(parests_all, [:code])
    CSV.write("results_firstdiff/select_country_siler_est_results.csv", parests_all)

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
@time map_dyn = optimize(log_siler_dyn_firstdiff(country_lm_data, country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
print("Estimated MAP for "*code)
# Estimate by MCMC
@time chain_dyn = sample(log_siler_dyn_firstdiff(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
    niters, nthreads, init_params = map_dyn.values.array)
print("Sampled posterior for "*code)
# Plot some example posterior distributions
display(chain_dyn)
plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]",
    "lc[1]", "ld[1]", "lσ[1]"]])
plot!(margin=8Plots.mm)
savefig("results_firstdiff/"*code*"/"*code*"_first_period_posteriors.pdf")
plot(chain_dyn[["σ_par[1]", "σ_par[2]", "σ_par[3]", "σ_par[4]", "σ_par[5]", "σ_par[6]"]])
plot!(margin=8Plots.mm)
savefig("results_firstdiff/"*code*"/"*code*"_rw_variance_posteriors.pdf")
plot(chain_dyn[["α_pars[1]", "α_pars[2]", "α_pars[3]", "α_pars[4]", "α_pars[5]", "α_pars[6]"]])
plot!(margin=8Plots.mm)
savefig("results_firstdiff/"*code*"/"*code*"_rw_drift_posteriors.pdf")
plot(chain_dyn[["β_pars[1]", "β_pars[2]", "β_pars[3]", "β_pars[4]", "β_pars[5]", "β_pars[6]"]])
plot!(margin=8Plots.mm)
savefig("results_firstdiff/"*code*"/"*code*"_rw_firstdiff_posteriors.pdf")

# Extract and plot results
parests_dyn = extract_variables(chain_dyn, country_years, log_pars = true,
    σ_pars = true, ext = true)
p1 = plot_siler_params(parests_dyn)
p_title = plot(title = "Dynamic (first diff) Siler parameters "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("results_firstdiff/"*code*"/"*code*"_param_estimates.pdf")
# Plot the time series parameters
p2 = plot_ts_params(parests_dyn)
p_title = plot(title = "Siler time series (first diff) parameters "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("results_firstdiff/"*code*"/"*code*"_ts_params.pdf")

plot_fit_year(parests_dyn, country_m_data[1], country_years[1], log_vals = false, col = 1)
plot_fit_year!(parests_dyn, country_m_data[T], country_years[T], log_vals = false, col = 2)
plot!(title = "Siler model fit "*string(code))
savefig("results_firstdiff/"*code*"/"*code*"_model_fit.pdf")

# Store and export results
CSV.write("results_firstdiff/"*code*"/"*code*"_est_results.csv", parests_dyn)
insertcols!(parests_dyn, 1, :code => repeat([Symbol(code)], nrow(parests_dyn)) )
insertcols!(parests_dyn, 2, :name => repeat([name], nrow(parests_dyn)) )

parests_all = vcat(parests_all, parests_dyn)
sort!(parests_all, [:code])
CSV.write("results_firstdiff/select_country_siler_est_results.csv", parests_all)


parests_dict[Symbol(code)] = parests_dyn
