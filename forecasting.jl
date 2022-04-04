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
nahead = 6 # Number of periods to forecast ahead

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

# Compute the model implied LE and H for the in-sample periods
df_post = compute_LE_post(df_post, years, nahead, spec = :Bergeron)



"""
Compute predictive distribution of forecasts
"""
# Compute df_pred which extends df_post with predictions
df_pred = compute_forecasts(df_post, nahead, ndraws, years; spec = :Colchero)
# Summarize past and forecast variables
past_years = years
fut_years = Int.(maximum(years) .+ 5.0.*(1:nahead))
parests_pred_col = extract_forecast_variables(df_pred, past_years, fut_years,
    log_pars = true, spec = :Bergeron, model_vers = :i2drift)
CSV.write(folder*"/siler_"*model*"_preds_col.csv", parests_pred_col)
# Plot Siler parameter forecasts
plot_siler_params(parests_pred_col, forecasts = true)
savefig(folder*"/siler_"*model*"_param_pred_col.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred_col, model_vers = :i2drift, forecasts = true)
savefig(folder*"/siler_"*model*"_ts_pred_col.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred_col, forecasts = true, bands = false)
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


"""
Check the life expectancy at later ages
"""
# Which extra ages do we want?
eval_ages = Int.(0:10:90)
# Add extra rows for these remianing life expectancies
rle_df = DataFrame(year = unique(parests_pred_col.year[parests_pred_col.year.>0]), LE0 = 0.)
for eval_age in eval_ages
    rle_df[:,Symbol("LE"*string(eval_age))] .= NaN
    for yy in rle_df.year
        year_df = parests_pred_col[parests_pred_col.year .== yy,:]
        params = SilerParam(b = year_df.median[year_df.parameter.==:b][1],
            B = year_df.median[year_df.parameter.==:B][1],
            c = year_df.median[year_df.parameter.==:c][1],
            C = year_df.median[year_df.parameter.==:C][1],
            d = year_df.median[year_df.parameter.==:d][1])

        rle_df[rle_df.year .== yy,Symbol("LE"*string(eval_age))] .=
            LE(params, eval_age, spec = :Colchero)
    end
end

plot(rle_df.year, rle_df.LE0, label = "Age 0",
    legend = :topleft, xlabel = "Year", ylabel = "Remaining LE")
plot!(rle_df.year, rle_df.LE10, label = "Age 10")
plot!(rle_df.year, rle_df.LE20, label = "Age 20")
plot!(rle_df.year, rle_df.LE30, label = "Age 30")
plot!(rle_df.year, rle_df.LE40, label = "Age 40")
plot!(rle_df.year, rle_df.LE50, label = "Age 50")
plot!(rle_df.year, rle_df.LE60, label = "Age 60")
plot!(rle_df.year, rle_df.LE70, label = "Age 70")
plot!(rle_df.year, rle_df.LE80, label = "Age 80")
plot!(rle_df.year, rle_df.LE90, label = "Age 90")
vline!([2020,2020], linestyle=:dot, color = :black, label = false)
plot!(size = (600,600))
