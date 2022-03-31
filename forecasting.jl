"""
Use MCMC from Turing.jl to estimate log Siler model on log mortality data
"""

if occursin("jashwin", pwd())
    cd("C://Users/jashwin/Documents/GitHub/MortalityEstimation/")
else
    cd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
end


# Import libraries.
using Turing, StatsPlots, Random, Optim, StatsBase, LinearAlgebra, Optim
using TruncatedDistributions, PDMats
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings, JLD2, ProgressMeter

include("src/MortalityEstimation.jl")
gr()

code = "Best Practice"
folder = "figures/benchmark/"
model = "i2"
ndraws = 10 # Number of draws to approximate each future shock
nahead = 4 # Number of periods to forecast ahead

# Import the full posterior
@load folder*"siler_"*model*"_chain.jld2" chain_i2
chain_post = deepcopy(chain_i2)
df_post = DataFrame(chain_post)
# Import the already extracted parameter estimates
parests_df_col =  CSV.read(folder*"siler_"*model*"_params_col.csv", DataFrame, ntasks = 1)
years = unique(parests_df_col.year[parests_df_col.year .> 0])
# Re-extract the parameters as a sense check
parests_df_col_check = extract_variables(chain_i2, years, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
plot_siler_params(parests_df_col_check)
plot_ts_params(parests_df_col_check, model_vers = :i2drift)



"""
 Compute the model implied LE and H for the in-sample periods
"""
# Add extra columns for model impled LE and H
LEs = Symbol.("LE[".*string.(1:length(years)+nahead).*"]")
Hs = Symbol.("H[".*string.(1:length(years)+nahead).*"]")
impl_vars = vcat(LEs, Hs)
df_post = hcat(df_post,DataFrame(NaN.*zeros(nrow(df_post), length(impl_vars)), impl_vars))
# Loop through each period
prog = Progress(length(years), desc = "Calculating model implied LE and H: ")
for ii in 1:length(years)
    params = particles2params(df_post, ii, log_pars = true)
    df_post[:, LEs[ii]] = LE.(params, [0.0], spec = :Colchero)
    df_post[:, Hs[ii]] = H.(params, [0.0], spec = :Colchero)
    next!(prog)
end




"""
Compute predictive distribution of forecasts
"""
# Compute df_pred which extends df_post with predictions
df_pred = compute_forecasts(df_post, nahead, ndraws, years; spec = :Colchero)
# Summarize past and forecast variables
past_years = years
fut_years = Int.(maximum(years) .+ 5.0.*(1:nahead))
parests_pred_col = extract_forecast_variables(df_pred, past_years, fut_years,
    log_pars = true, spec = :Colchero, model_vers = :i2drift)
CSV.write(folder*"/siler_"*model*"_preds_col.csv", parests_pred_col)
# Plot Siler parameter forecasts
plot_siler_params(parests_pred_col, forecasts = true)
savefig(folder*"/siler_"*model*"_param_pred_col.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred_col, model_vers = :i2drift, forecasts = true)
savefig(folder*"/siler_"*model*"_ts_pred_col.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred_col, forecasts = true)
savefig(folder*"/siler_"*model*"_leh_pred_col.pdf")
# Forecast decomposition of future LE and H
decomp_pred_col = create_decomp(parests_pred_col[parests_pred_col.year .> 1985,:],
    spec = :Colchero, eval_age = 0)
le_p = plot_decomp(decomp_pred_col, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_pred_col, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Future decomposition Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig(folder*"/siler_"*model*"_decomp_pred_col.pdf")
