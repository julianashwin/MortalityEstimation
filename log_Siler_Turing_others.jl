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
    "JPN", "KOR", "NLD", "NZL_NM", "NOR", "PRT", "USA"]
country_codes = ["AUT", "BGR", "BLR", "CHL", "CZE", "DEU", "DEUTW", "DNK", "EST", "HRV",
    "HUN", "IRL", "ISR", "LTU", "LUX", "LVA", "POL", "RUS", "SVK", "SVN", "TWN", "UKR"]

country_codes = ["AUS", "CAN", "CHE", "BEL", "ESP", "FIN", "FRA", "GBR", "GRC", "HKG", "ITA", "ISL",
    "JPN", "KOR", "NLD", "NZL_NM", "NOR", "PRT", "USA", "AUT", "BGR", "BLR", "CHL", "CZE", "DEU",
    "DEUTW", "DNK", "EST", "HRV","HUN", "IRL", "ISR", "LTU", "LUX", "LVA", "POL", "RUS", "SVK",
    "SVN", "TWN", "UKR"]


#unique(mort_df.code)[.!(in.(unique(mort_df.code),[country_codes]))]

select_df = mort_df[in.(mort_df.code, [country_codes]),:]

## Set model and folder to save results
folder = "countries"
model = "i2"



"""
Check data looks sensible for a single country
"""
## Data prep for single coujntry
code = country_codes[18]
#country_df = mort_df[(mort_df.code .== code), :]
country_df = select_df[select_df.code .== code,:]
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
Also save the decompositions
"""
# Loop through each country to compute the predictions and decompose
ndraws = 10 # Number of draws to approximate each future shock
nahead = 6 # Number of periods to forecast ahead
for code in country_codes
    print("Save decomposition for "*code*"\n")
    country_df = select_df[(select_df.code .== code), :]
    country_years = unique(country_df.year)

    @load "figures/"*folder*"/"*code*"_"*model*"_chain.jld2" chain_i2

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
    # Extract decomp
    decomp_df = create_decomp(parests_pred; spec = :Bergeron, eval_age = 0)
    CSV.write("figures/"*folder*"/"*code*"_"*model*"_decomp_pred.csv", decomp_df)
    # Extract gradients
    LEgrad_df = compute_LEgrad_df(decomp_df; ages = 0:140)
    CSV.write("figures/"*folder*"/"*code*"_"*model*"_LEgrads.csv", LEgrad_df)


end





"""
Out-of-sample forecasts
"""

folder = "countries/held-out"
model = "i2"

for code in country_codes[[11,12,13]]
    print("Working on oos forecasts for "*code)
    # Extract and convert relevant data into correct form
    country_df = select_df[(select_df.code .== code), :]
    country_df = country_df[country_df.year .>=1900,:]
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
    if T >= 18
        for held_out = 2:2:10
            # Adjust if you don't want every period
            periods = Int.(1:T-held_out)
            years_selected = Int.(round.(country_years[periods]))

            ## Find some starting points
            # MAP estimate for multiple independent models on log mortality
            #map_i2 = optimize(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), MAP(), LBFGS(),
            #    Optim.Options(iterations=60_000, allow_f_increases=true))
            #map_i2_vals =  map_i2.values.array
            # Alternatively, we can simulate from the prior and start there
            prior_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), Prior(), 5000)
            df_prior = DataFrame(prior_i2)
            insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
            showtable(df_prior)
            prior_i2_vals = df_prior[1,3:(end-1)]

            ## Estimate the model
            niters = 1000
            nchains = 4
            # MCMC sampling
            chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
                niters, nchains, init_params = prior_i2_vals)
            display(chain_i2)
            @save "figures/"*folder*"/"*code*"_siler_"*model*"_chain_"*string(years_selected[T-held_out])*".jld2" chain_i2

        end
    end
end
