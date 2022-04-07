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
parests_df =  CSV.read(folder*"siler_"*model*"_params_ber.csv", DataFrame, ntasks = 1)
years = unique(parests_df.year[parests_df.year .> 0])
# Re-extract the parameters as a sense check
parests_df_check = extract_variables(chain_post, years, log_pars = true,
    model_vers = :i2drift, spec = :Bergreon)
plot_siler_params(parests_df_check)
plot_ts_params(parests_df_check, model_vers = :i2drift)

# Compute the model implied LE and H for the in-sample periods
df_post = compute_LE_post(df_post, years, nahead, spec = :Bergeron)


# Import the empirical data
mort_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1)
bp_df = mort_df[mort_df.best_practice.==1,:]
sort!(bp_df, [:year, :age])



"""
Casual little decomposition
"""
# Decompose historical chandes to make sure the results look right
decomp_df = create_decomp(parests_df; spec = :Bergeron, eval_age = 0)
le_p = plot_decomp(decomp_df, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px,
    yticks = false, xticks = false)
plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400))
display(plot!(left_margin = 10Plots.px, bottom_margin = 10Plots.px))
savefig(folder*"/siler_"*model*"_decomp_ber.pdf")

decomp_post60 = decomp_df[decomp_df.year .>= 1960,:]
sum(decomp_post60.ΔC.*decomp_post60.LE_C)/sum(decomp_post60.ΔLE_mod)


# Decompose with Colchero specification for comparison
parests_df_col = extract_variables(chain_post, years, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
decomp_df_col = create_decomp(parests_df_col; spec = :Colchero, eval_age = 0)
le_p = plot_decomp(decomp_df_col, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_col, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Colcher parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px,
    yticks = false, xticks = false)
plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400))
display(plot!(left_margin = 10Plots.px, bottom_margin = 10Plots.px))
savefig(folder*"/siler_"*model*"_decomp_col.pdf")



"""
Compute predictive distribution of forecasts
"""
# Compute df_pred which extends df_post with predictions
df_pred = compute_forecasts(df_post, nahead, ndraws, years; spec = :Bergeron)
# Summarize past and forecast variables
past_years = years
fut_years = Int.(maximum(years) .+ 5.0.*(1:nahead))
parests_pred = extract_forecast_variables(df_pred, past_years, fut_years,
    log_pars = true, spec = :Bergeron, model_vers = :i2drift)
#CSV.write(folder*"/siler_"*model*"_preds_col.csv", parests_pred_col)
# Plot Siler parameter forecasts
plot_siler_params(parests_pred, forecasts = true)
#savefig(folder*"/siler_"*model*"_param_pred_col.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred, model_vers = :i2drift, forecasts = true)
#savefig(folder*"/siler_"*model*"_ts_pred_col.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred, forecasts = true, bands = true)
#savefig(folder*"/siler_"*model*"_leh_pred_col.pdf")
# Forecast decomposition of future LE and H
decomp_pred = create_decomp(parests_pred[parests_pred_col.year .> 1985,:],
    spec = :Bergeron, eval_age = 0)
le_p = plot_decomp(decomp_pred, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
vline!([2020.5,2020.5], linestyle=:dot, color = :red, label = false)
h_p = plot_decomp(decomp_pred, :H)
h_p = plot!(title = "Lifespan Inequality")
vline!([2020.5,2020.5], linestyle=:dot, color = :red, label = false)
p_title = plot(title = "Future decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400))
display(plot!(left_margin = 10Plots.px, bottom_margin = 10Plots.px))
savefig(folder*"/siler_"*model*"_decomp_pred.pdf")


decomp_post2020 = decomp_pred[decomp_pred.year .>= 2020,:]
sum(decomp_post2020.ΔC.*decomp_post2020.LE_C)/sum(decomp_post2020.ΔLE_mod)


"""
Check the life expectancy at later ages
"""
# Which extra ages do we want?
eval_ages = Int.(0:10:90)
# Add extra rows for these remianing life expectancies
rle_df = DataFrame(year = unique(parests_pred_col.year[parests_pred_col.year.>0]), LE0 = 0.)
for eval_age in eval_ages
    rle_df[:,Symbol("LE"*string(eval_age))] .= NaN
    rle_df[:,Symbol("ex"*string(eval_age))] .= NaN
    for yy in rle_df.year
        year_df = parests_pred_col[parests_pred_col.year .== yy,:]
        params = SilerParam(b = year_df.median[year_df.parameter.==:b][1],
            B = year_df.median[year_df.parameter.==:B][1],
            c = year_df.median[year_df.parameter.==:c][1],
            C = year_df.median[year_df.parameter.==:C][1],
            d = year_df.median[year_df.parameter.==:d][1])

        rle_df[rle_df.year .== yy,Symbol("LE"*string(eval_age))] .=
            LE(params, eval_age, spec = :Bergeron)
        if yy in unique(bp_df.year)
            rle_df[rle_df.year .== yy,Symbol("ex"*string(eval_age))] .=
                bp_df[(bp_df.year .== yy) .& (bp_df.age .== eval_age),:ex_f]
        end
    end
end

plot(rle_df.year, rle_df.LE0, label = "Age 0", color = 1,
    legend = :outerright, xlabel = "Year", ylabel = "Remaining LE")
scatter!(rle_df.year, rle_df.ex0, label = false, markershape = :cross, color = 1)
plot!(rle_df.year, rle_df.LE10, label = "Age 10", color = 2)
scatter!(rle_df.year, rle_df.ex10, label = false, markershape = :cross, color = 2)
plot!(rle_df.year, rle_df.LE20, label = "Age 20", color = 3)
scatter!(rle_df.year, rle_df.ex20, label = false, markershape = :cross, color = 3)
plot!(rle_df.year, rle_df.LE30, label = "Age 30", color = 4)
scatter!(rle_df.year, rle_df.ex30, label = false, markershape = :cross, color = 4)
plot!(rle_df.year, rle_df.LE40, label = "Age 40", color = 5)
scatter!(rle_df.year, rle_df.ex40, label = false, markershape = :cross, color = 5)
plot!(rle_df.year, rle_df.LE50, label = "Age 50", color = 6)
scatter!(rle_df.year, rle_df.ex50, label = false, markershape = :cross, color = 6)
plot!(rle_df.year, rle_df.LE60, label = "Age 60", color = 7)
scatter!(rle_df.year, rle_df.ex60, label = false, markershape = :cross, color = 7)
plot!(rle_df.year, rle_df.LE70, label = "Age 70", color = 8)
scatter!(rle_df.year, rle_df.ex70, label = false, markershape = :cross, color = 8)
plot!(rle_df.year, rle_df.LE80, label = "Age 80", color = 9)
scatter!(rle_df.year, rle_df.ex80, label = false, markershape = :cross, color = 9)
plot!(rle_df.year, rle_df.LE90, label = "Age 90", color = 10)
scatter!(rle_df.year, rle_df.ex90, label = false, markershape = :cross, color = 10)
vline!([2020,2020], linestyle=:dot, color = :black, label = false, ylims = (0,100))
plot!(size = (600,300))
savefig(folder*"/siler_"*model*"_rle_dist_pred.pdf")




"""
Plot the gradient of each parameter over time - are they becoming more effective?
"""
function plot_LEgrads(decomp_df)
    le_p = plot(layout = (2,2))
    le_p = plot!(decomp_df.year, decomp_df.LE_b, label = "b", subplot = 1)
    le_p = hline!([0,0], color = :black, linestyle = :dash, label = false, subplot = 1)
    le_p = plot!(decomp_df.year, decomp_df.LE_B, label = "B", subplot = 2)
    le_p = hline!([0,0], color = :black, linestyle = :dash, label = false, subplot = 2)
    le_p = plot!(decomp_df.year, decomp_df.LE_c, label = "c", subplot = 3)
    le_p = hline!([0,0], color = :black, linestyle = :dash, label = false, subplot = 3)
    le_p = plot!(decomp_df.year, decomp_df.LE_C, label = "C", subplot = 4)
    le_p = hline!([0,0], color = :black, linestyle = :dash, label = false, subplot = 4)
    #le_p = plot!(decomp_df.year, decomp_df.LE_d, label = "d", subplot = 1)

    return le_p
end

decomp_pred = create_decomp(parests_pred, spec = :Bergeron, eval_age = 0)
pre_2020 = decomp_pred.year .<= 2018
post_2020 = decomp_pred.year .>= 2018
# Bergeron
bB_plt = plot(layout = (2,2), legend = false)
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_b[pre_2020], title = L"LE_b", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_b[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,7))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_B[pre_2020], title = L"LE_B", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_B[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,3))
plot!(decomp_pred.year[pre_2020], decomp_pred.H_b[pre_2020], title = L"H_b", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.H_b[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, ylims = (-0.1,0), xmirror = true)
plot!(decomp_pred.year[pre_2020], decomp_pred.H_B[pre_2020], title = L"H_B", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.H_B[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, ylims = (-0.05,0), xmirror = true)
plot!(size = (600,300))
savefig(folder*"siler_"*model*"_LEgrad_bB.pdf")


cC_plt = plot(layout = (2,2), legend = false)
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_c[pre_2020], title = L"LE_c", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_c[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,75))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_C[pre_2020], title = L"LE_C", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_C[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,1))
plot!(decomp_pred.year[pre_2020], decomp_pred.H_c[pre_2020], title = L"H_c", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.H_c[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, ylims = (-2.4,0), xmirror = true)
plot!(decomp_pred.year[pre_2020], decomp_pred.H_C[pre_2020], title = L"H_C", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.H_C[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, ylims = (-0.002,0), xmirror = true)
plot!(size = (600,300))
savefig(folder*"siler_"*model*"_LEgrad_cC.pdf")


# What if we set d to zero?
decomp_today = decomp_df[decomp_df.year .== 2018,[:B,:b,:C,:c,:d,:σ, :LE_mod]]
param_today = SilerParam(b= decomp_today.b[1], B = decomp_today.B[1], c = decomp_today.c[1],
    C = decomp_today.C[1], d = decomp_today.d[1])
param_dzero = SilerParam(b= decomp_today.b[1], B = decomp_today.B[1], c = decomp_today.c[1],
    C = decomp_today.C[1], d = 0.)

LE(param_today, 0.0, spec = :Bergeron)
LE(param_dzero, 0.0, spec = :Bergeron)
