

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
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings, JLD2


include("src/MortalityEstimation.jl")

# Import and extract the best practice series
all_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1)
bp_df = all_df[(all_df.best_practice .== 1) .& (all_df.year .>= 1900), :]
sort!(bp_df, [:year, :age])


code = "Best Practice"
folder = "benchmark"
model = "i2"

@load "figures/"*folder*"/siler_i2_chain.jld2" chain_i2

years_selected = unique(bp_df.year)

ndraws = 10 # Number of draws to approximate each future shock
nahead = 6 # Number of periods to forecast ahead
df_post = DataFrame(chain_i2)
# Compute the model implied LE and H for the in-sample periods
df_post = compute_LE_post(df_post, years_selected, nahead, spec = :Bergeron, model_vers = :i2drift)
# Compute df_pred which extends df_post with predictions
df_pred = compute_forecasts(df_post, nahead, ndraws, years_selected, spec = :Bergeron, model_vers = :i2drift)
# Summarize past and forecast variables
past_years = years_selected
fut_years = Int.(maximum(years_selected) .+ 5.0.*(1:nahead))
parests_pred = extract_forecast_variables(df_pred, past_years, fut_years,
    log_pars = true, spec = :Bergeron, model_vers = :i2drift)

# Plot Siler parameter forecasts
plot_siler_params(parests_pred, forecasts = true)
savefig("figures/"*folder*"/siler_"*model*"_param_pred.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred, model_vers = :i2drift, forecasts = true)
savefig("figures/"*folder*"/siler_"*model*"_ts_pred.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred, forecasts = true, bands = true)
savefig("figures/"*folder*"/siler_"*model*"_leh_pred.pdf")
# Plot forecast decomposition
decomp_pred = create_decomp(parests_pred[parests_pred.year .> 1900,:],
    spec = :Bergeron, eval_age = 0, forecasts = true)
LE_p = plot_decomp(decomp_pred, :LE, forecasts = true)
LE_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_pred, :h, forecasts = true)
h_p = plot!(title = "Lifespan Equality")
H_p = plot_decomp(decomp_pred, :H, forecasts = true)
H_p = plot!(legend = false, title = "Lifespan Inequality")
p_title = plot(title = "Future decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
plot(p_title, LE_p, h_p, H_p, layout = @layout([A{0.01h}; B C D]), size = (1000,400))
display(plot!(left_margin = 20Plots.px, bottom_margin = 15Plots.px))
savefig("figures/"*folder*"/siler_"*model*"_decomp_pred.pdf")
