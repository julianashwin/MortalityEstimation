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
all_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1)

# Best practice data
mort_df = deepcopy(all_df)  # [(all_df.year .>= 1700), :]
sort!(mort_df, [:code, :year, :age])
bp_df = mort_df[(mort_df.best_practice .== 1), :]
plot(bp_df.age, bp_df.mx, group = bp_df.year, legend = :top)
sort!(bp_df, [:year, :age])

# Restrict to G7 post 1900 for now
country_codes = ["AUS", "CAN", "CHE", "BEL", "ESP", "FIN", "FRA", "GBR", "GRC", "HKG", "ITA", "ISL",
    "JPN", "KOR", "NZL_NM", "NOR", "PRT", "USA"]
select_df = mort_df[in.(mort_df.code, [country_codes]),:]

## Set model and folder to save results
folder = "countries"
model = "i2"


"""
Check data looks sensible for a single country
"""
## Data prep for single coujntry
code = "KOR"
#country_df = mort_df[(mort_df.code .== code), :]
country_df = select_df[select_df.code .== "KOR",:]
# Need to remove any zeros
country_df[country_df[:,:mx_f] .== 0.0,:mx_f] .=  minimum(country_df[country_df[:,:mx_f] .> 0.0,:mx_f])
# Check data looks sensible
plot(country_df.age, country_df.mx_f, group = country_df.year, legend = :top)
# Convert this into a matrix of mortality rates over time, age and year vectors
country_m_data = chunk(country_df.mx_f, 110)
plot(country_m_data, legend = false)
country_lm_data = [log.(m_dist) for m_dist in country_m_data]
plot(country_lm_data, legend = false)
country_ages = Int64.(0:maximum(country_df.age))
country_years = unique(country_df.year)
T = length(country_lm_data)

@assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
@assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"




"""
Estimate the I(2) model on each of the selected countries
"""
# Number of iterations and chains for sampler
niters = 1250
nchains = 4
ndraws = 10 # Number of draws to approximate each future shock
nahead = 6 # Number of periods to forecast ahead
# Loop through the selected codes
for code in country_codes
    print("Working on model for "*code)
    # Extract and convert relevant data into correct form
    country_df = select_df[(select_df.code .== code), :]
    sort!(select_df, [:year, :age])
    # Suprisingly, we actually have some zeros here for small countries (e.g. ISL)
    country_df.mx_f[country_df.mx_f .== 0.0] .=  minimum(country_df.mx_f[country_df.mx_f .> 0.0])
    # Get data into right format
    country_m_data = chunk(country_df.mx_f, 110)
    country_lm_data = [log.(m_dist) for m_dist in country_m_data]
    country_ages = Int64.(0:maximum(country_df.age))
    country_years = unique(country_df.year)
    T = length(country_lm_data)
    name = country_df.name[1]
    # Simulate from prior as starting point
    prior_i2 = sample(log_siler_dyn_i2drift(country_lm_data, country_ages), Prior(), 5000)
    df_prior = DataFrame(prior_i2)
    insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
    prior_i2_vals = df_prior[1,3:(end-1)]
    # Estimate posterior with NUTS sampler
    chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
        niters, nchains, init_params = prior_i2_vals)
    # Save chain
    @save "figures/"*folder*"/"*code*"_"*model*"_chain.jld2" chain_i2
    # Extract some summary statistics
    parests_i2 = extract_variables(chain_i2, country_years, log_pars = true,
        model_vers = :i2drift, spec = :Bergeron)
    CSV.write("figures/"*folder*"/"*code*"_"*model*"_params.csv", parests_i2)
    # Plot the parameters
    display(plot_siler_params(parests_i2))
    savefig("figures/"*folder*"/"*code*"_"*model*"_params.pdf")
    plot_ts_params(parests_i2, model_vers = :i2drift)
    savefig("figures/"*folder*"/"*code*"_"*model*"_ts_params.pdf")
    # Compute decomposition and plot
    decomp_df = create_decomp(parests_i2; spec = :Bergeron, eval_age = 0)
    le_p = plot_decomp(decomp_df, :LE)
    le_p = plot!(legend = false, title = "Life Expectancy")
    h_p = plot_decomp(decomp_df, :H)
    h_p = plot!(title = "Lifespan Inequality")
    p_title = plot(title = "Historical decomposition Siler Bergeron parameters "*string(code),
        grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
    plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400))
    display(plot!(left_margin = 10Plots.px, bottom_margin = 10Plots.px))
    savefig("figures/"*folder*"/"*code*"_"*model*"_decomp.pdf")
    # Compute some forecasts
    df_post = DataFrame(chain_i2)
    # Compute the model implied LE and H for the in-sample periods
    df_post = compute_LE_post(df_post, country_years, nahead, spec = :Bergeron)
    # Compute df_pred which extends df_post with predictions
    df_pred = compute_forecasts(df_post, nahead, ndraws, country_years, spec = :Bergeron)
    # Extract forecasts
    parests_pred = extract_forecast_variables(df_pred, country_years, Int.(maximum(country_years) .+ 5.0.*(1:nahead)),
        log_pars = true, spec = :Bergeron, model_vers = :i2drift)
    CSV.write("figures/"*folder*"/"*code*"_"*model*"_preds.csv", parests_pred)
    # Plot forecasts
    plot_siler_params(parests_pred, forecasts = true)
    savefig("figures/"*folder*"/"*code*"_"*model*"_param_pred.pdf")
    plot_ts_params(parests_pred, model_vers = :i2drift, forecasts = true)
    savefig("figures/"*folder*"/"*code*"_"*model*"_ts_pred.pdf")
    plot_LE_H(parests_pred, forecasts = true, bands = true)
    savefig("figures/"*folder*"/"*code*"_"*model*"_leh_pred.pdf")

end






















"""
Preliminary checks and illustrative estimation of first and last period
"""
## Data prep for single coujntry
country_df = select_df[(select_df.code .== "UKR"), :]
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
parests_indep = extract_variables(chain_indep, years_selected, log_pars = false,
    σ_pars = false)
plot_siler_params(parests_indep)
# Visualise model fit
plot_fit_year(parests_indep, country_m_data[1], years_selected[1], log_vals = false, col =1)
plot_fit_year!(parests_indep, country_m_data[T], years_selected[end], log_vals = false, col = 2)



"""
Estimate Siler model for each country in G7
"""
## Set some estimation options
niters = 1000
nthreads = 4

parests_all = DataFrame(code = Symbol[], name = String[], parameter = Symbol[], year = Int64[],
    mean = Float64[], min = Float64[], median = Float64[], max = Float64[],
    nmissing  = Int64[], eltype = DataType[], std = Float64[],
    pc975 = Float64[], pc025 = Float64[], pc85 = Float64[], pc15 = Float64[],
    pc75 = Float64[], pc25 = Float64[], mcse = Float64[], ess = Float64[], rhat = Float64[],
    ess_per_sec = Float64[])

parests_dict = Dict{Symbol, DataFrame}()


for code in unique(G7_df.code)[11:36]
    print("Working on model for "*code)
    # Extract and convert relevant data into correct form
    country_df = G7_df[(G7_df.code .== code), :]

    # Suprisingly, we actually have some zeros here for small countries (e.g. ISL)
    country_df.mx[country_df.mx .== 0.0] .=  minimum(country_df.mx[country_df.mx .> 0.0])

    # Get data into right format
    country_m_data = chunk(country_df.mx, 110)
    country_lm_data = [log.(m_dist) for m_dist in country_m_data]
    country_ages = Int64.(0:maximum(country_df.age))
    country_years = unique(country_df.year)
    T = length(country_lm_data)
    name = country_df.name[1]

    @assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
    @assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"

    # Sample from prior as a sense check
    chain_prior = sample(log_siler_dyn_ext(country_lm_data, country_ages), Prior(), 5000)
    df_prior = DataFrame(chain_prior)
    insert!.(eachcol(df_prior), 1, vcat([0,0],mean.(eachcol(df_prior[:,3:end]))))
    prior_means = mean.(eachcol(df_prior[:,3:end-1]))

    # MAP estimate to initialise MCMC
    @time map_dyn = optimize(log_siler_dyn_ext(country_lm_data, country_ages), MAP(), LBFGS(),
        Optim.Options(iterations=50_000, allow_f_increases=true))
    print("Estimated MAP for "*code)
    # Estimate by MCMC
    @time chain_dyn = sample(log_siler_dyn_ext(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
        niters, nthreads, init_params = map_dyn.values.array)
    print("Sampled posterior for "*code)
    # Plot some example posterior distributions
    display(chain_dyn)
    plot(chain_dyn[["lB[1]", "lb[1]", "lC[1]",
        "lc[1]", "ld[1]", "lσ[1]"]])
    plot!(margin=8Plots.mm)
    savefig("results/"*code*"/"*code*"_first_period_posteriors.pdf")
    plot(chain_dyn[["σ_par[1]", "σ_par[2]", "σ_par[3]", "σ_par[4]", "σ_par[5]", "σ_par[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results/"*code*"/"*code*"_rw_variance_posteriors.pdf")
    plot(chain_dyn[["α_pars[1]", "α_pars[2]", "α_pars[3]", "α_pars[4]", "α_pars[5]", "α_pars[6]"]])
    plot!(margin=8Plots.mm)
    savefig("results/"*code*"/"*code*"_rw_variance_posteriors.pdf")
    # Extract and plot results
    parests_dyn = extract_variables(chain_dyn, country_years, log_pars = true,
        σ_pars = true, ext = true)
    p1 = plot_siler_params(parests_dyn)
    p_title = plot(title = "Dynamic Siler parameters "*string(code), grid = false, showaxis = false,
        bottom_margin = -10Plots.px, yticks = false, xticks = false)
    display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
    savefig("results/"*code*"/"*code*"_param_estimates.pdf")
    # Plot the time series parameters
    p2 = plot_ts_params(parests_dyn)
    p_title = plot(title = "Siler time series parameters "*string(code), grid = false, showaxis = false,
        bottom_margin = -10Plots.px, yticks = false, xticks = false)
    display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
    savefig("results/"*code*"/"*code*"_ts_params.pdf")
    # Plot fit in selected years
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
    CSV.write("results/other_country_siler_est_results.csv", parests_all)

    parests_dict[Symbol(code)] = parests_dyn

end
