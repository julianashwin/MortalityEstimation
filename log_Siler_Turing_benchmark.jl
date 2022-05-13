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
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings, JLD2


include("src/MortalityEstimation.jl")


## Import the mortality data
all_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1)

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
Set up a single country/case to run through examples
"""
## Data prep for single coujntry
code = "Best Practice"
#code = "KOR"
folder = "benchmark"
#country_df = mort_df[(mort_df.code .== code), :]
country_df = bp_df
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


lifespan_df = DataFrame(year = unique(country_df.year), Lstar = 0.)
for yy in 1:nrow(lifespan_df)
    year = lifespan_df.year[yy]
    year_df = country_df[country_df.year .== year,:]
    plot(year_df.age, year_df.lx_f)
    Lstar_emp = year_df.age[year_df.lx_f .<= 0.001]
    if length(Lstar_emp) > 0
        Lstar_emp = minimum(Lstar_emp)
    else
        Lstar_emp = 110
    end
    lifespan_df.Lstar[yy] = Lstar_emp
end

"""
Multiple independent Siler models
"""
# Define model to save results
model = "indep"
# Adjust if you don't want every period
periods = Int.(1:T)
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
niters = 1000
nchains = 4
# MCMC sampling
chain_indep = sample(log_siler_indep(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_indep_vals)
display(chain_indep)
@save "figures/"*folder*"/siler_"*model*"_chain.jld2" chain_indep

## Plot Siler parameters
# Display the parameters with Bergeron specification
parests_indep_ber = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Bergeron)
p1 = plot_siler_params(parests_indep_ber)
p_title = plot(title = "Multiple independent Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_ber.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_ber.csv", parests_indep_ber)
# With Colchero specification
parests_indep_col = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Colchero)
p1 = plot_siler_params(parests_indep_col)
p_title = plot(title = "Multiple independent Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_col.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_col.csv", parests_indep_col)
# With Scott specification
parests_indep_sco = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Scott)
p1 = plot_siler_params(parests_indep_sco)
p_title = plot(title = "Multiple independent Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_sco.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_sco.csv", parests_indep_sco)

## Plot decomposition
# Bergeron specification
decomp_df_ber = create_decomp(parests_indep_ber; spec = :Bergeron, eval_age = 0)
le_p = plot_decomp(decomp_df_ber, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_ber, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig("figures/"*folder*"/siler_"*model*"_decomp_ber.pdf")
# Colchero specification
decomp_df_col = create_decomp(parests_indep_col; spec = :Colchero, eval_age = 0)
le_p = plot_decomp(decomp_df_col, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_col, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig("figures/"*folder*"/siler_"*model*"_decomp_col.pdf")
# Scott specification
decomp_df_sco = create_decomp(parests_indep_sco; spec = :Scott, eval_age = 0)
le_p = plot_decomp(decomp_df_sco, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_sco, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig("figures/"*folder*"/siler_"*model*"_decomp_sco.pdf")


"""
Dynamic Siler model with parameters i(2) random walk with drift
"""
# Define model to save results
model = "i2drift"
folder = "benchmark"

## Simulate from the prior to start there
prior_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_i2)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
showtable(df_prior)
prior_i2_vals = df_prior[1,3:(end-1)]

## Estimate the model
niters = 1250
nchains = 4
# MCMC sampling
chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data, country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_i2_vals)
display(chain_i2)
@save "figures/"*folder*"/siler_i2_chain.jld2" chain_i2
#@load "figures/"*folder*"/siler_i2_chain.jld2" chain_i2

"""
Plot Siler parameters
"""
# Display the parameters with Bergeron specification
parests_i2_ber = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Bergeron)
p1 = plot_siler_params(parests_i2_ber)
p_title = plot(title = "I(2) Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_ber.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_ber.csv", parests_i2_ber)
# with Colchero specification
parests_i2_col = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
p1 = plot_siler_params(parests_i2_col)
p_title = plot(title = "I(2) Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_col.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_col.csv", parests_i2_col)


## Plot time series parameters
p2 = plot_ts_params(parests_i2_ber, model_vers = :i2drift)
p_title = plot(title = "I(2) Siler Bergeron ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_ts_params.pdf")


"""
Plot decomposition
"""
# Bergeron specification
decomp_df_ber = create_decomp(parests_i2_ber; spec = :Bergeron, eval_age = 0)
LE_p = plot_decomp(decomp_df_ber, :LE)
LE_p = plot!(legend = false, title = "Life Expectancy")
H_p = plot_decomp(decomp_df_ber, :H)
H_p = plot!(legend = false, title = "Lifespan Inequality")
h_p = plot_decomp(decomp_df_ber, :h)
h_p = plot!(title = "Lifespan Equality")
p_title = plot(title = "Historical decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
plot(p_title, LE_p, h_p, H_p, layout = @layout([A{0.01h}; B C D]), size = (1000,400))
display(plot!(left_margin = 15Plots.px, bottom_margin = 15Plots.px))
savefig("figures/"*folder*"/siler_"*model*"_decomp_ber.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_decomp_ber.csv", decomp_df_ber)

# Colchero specification
decomp_df_col = create_decomp(parests_i2_col; spec = :Colchero, eval_age = 0)
LE_p = plot_decomp(decomp_df_col, :LE)
LE_p = plot!(legend = false, title = "Life Expectancy")
H_p = plot_decomp(decomp_df_col, :H)
H_p = plot!(legend = false,title = "Lifespan Inequality")
h_p = plot_decomp(decomp_df_col, :h)
h_p = plot!(title = "Lifespan Equality")
p_title = plot(title = "Historical decomposition Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
plot(p_title, LE_p, h_p, H_p, layout = @layout([A{0.01h}; B C D]), size = (1000,400))
display(plot!(left_margin = 15Plots.px, bottom_margin = 15Plots.px))
savefig("figures/"*folder*"/siler_"*model*"_decomp_col.pdf")


"""
Forecasts
"""
# Options
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
CSV.write("figures/"*folder*"/siler_"*model*"_preds.csv", parests_pred)

# Plot Siler parameter forecasts
plot_siler_params(parests_pred, forecasts = true, bands = true)
savefig("figures/"*folder*"/siler_"*model*"_param_pred.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred, model_vers = :i2drift, forecasts = true)
savefig("figures/"*folder*"/siler_"*model*"_ts_pred.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred, forecasts = true, bands = true)
savefig("figures/"*folder*"/siler_"*model*"_leh_pred.pdf")
# Plot forecasts of model implied LE,Lstar and Lmed
plot_Ls(parests_pred, forecasts = true, bands = true)
savefig("figures/"*folder*"/siler_"*model*"_ls_pred.pdf")
scatter!(lifespan_df.year, lifespan_df.Lstar, label = false, markershape = :cross,
    subplot = 2, color = :green)

# Forecast decomposition of future LE and H
decomp_pred = create_decomp(parests_pred[parests_pred.year .> 1985,:],
    spec = :Bergeron, eval_age = 0, forecasts = true)
LE_p = plot_decomp(decomp_pred, :LE, forecasts = true)
LE_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_pred, :h, forecasts = true)
h_p = plot!(title = "Lifespan Equality")
Lstar_p = plot_decomp(decomp_pred, :Lstar, forecasts = true)
Lstar_p = plot!(legend = false, title = "Lifespan")
p_title = plot(title = "Future decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
plot(p_title, LE_p, Lstar_p, h_p, layout = @layout([A{0.01h}; B C D]), size = (1000,400))
display(plot!(left_margin = 20Plots.px, bottom_margin = 15Plots.px))
savefig("figures/"*folder*"/siler_"*model*"_decomp_pred.pdf")



"""
LE gradients over time
"""
## Prediction decomposition
decomp_pred = create_decomp(parests_pred, spec = :Bergeron, eval_age = 0, forecasts = true)
pre_2020 = decomp_pred.year .<= 2018
post_2020 = decomp_pred.year .>= 2018
## Infant
bB_plt = plot(layout = (3,2), legend = false)
# LE
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_b[pre_2020], title = L"LE_b", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_b[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,maximum(decomp_pred.LE_b)))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_B[pre_2020], title = L"LE_B", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_B[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,maximum(decomp_pred.LE_B)))
# Lstar
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_b[pre_2020], title = L"L^*_b", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_b[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, ylims = (0, maximum(decomp_pred.Lstar_b)))
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_B[pre_2020], title = L"L^*_B", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_B[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, ylims = (0, maximum(decomp_pred.Lstar_B)))
# h
plot!(decomp_pred.year[pre_2020], decomp_pred.h_b[pre_2020], title = L"h_b", subplot = 5)
plot!(decomp_pred.year[post_2020], decomp_pred.h_b[post_2020], linestyle = :dash,subplot = 5)
hline!([0,0], color = :black, linestyle = :solid, subplot = 5, ylims = (0, maximum(decomp_pred.h_b)))
plot!(decomp_pred.year[pre_2020], decomp_pred.h_B[pre_2020], title = L"h_B", subplot = 6)
plot!(decomp_pred.year[post_2020], decomp_pred.h_B[post_2020], linestyle = :dash,subplot = 6)
hline!([0,0], color = :black, linestyle = :solid, subplot = 6, ylims = (0, maximum(decomp_pred.h_B)))
plot!(size = (600,300))
savefig("figures/"*folder*"siler_"*model*"_LEgrad_bB.pdf")

## Senescent
cC_plt = plot(layout = (3,2), legend = false)
# LE
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_c[pre_2020], title = L"LE_c", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_c[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,maximum(decomp_pred.LE_c)))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_C[pre_2020], title = L"LE_C", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_C[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,maximum(decomp_pred.LE_C)))
# Lstar
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_c[pre_2020], title = L"L^*_c", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_c[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, ylims = (minimum(decomp_pred.Lstar_c), 0))
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_C[pre_2020], title = L"L^*_C", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_C[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, ylims = (0, maximum(decomp_pred.Lstar_C)))
# h
plot!(decomp_pred.year[pre_2020], decomp_pred.h_c[pre_2020], title = L"h_c", subplot = 5)
plot!(decomp_pred.year[post_2020], decomp_pred.h_c[post_2020], linestyle = :dash,subplot = 5)
hline!([0,0], color = :black, linestyle = :solid, subplot = 5, ylims = (0, maximum(decomp_pred.h_c)))
plot!(decomp_pred.year[pre_2020], decomp_pred.h_C[pre_2020], title = L"h_C", subplot = 6)
plot!(decomp_pred.year[post_2020], decomp_pred.h_C[post_2020], linestyle = :dash,subplot = 6)
hline!([0,0], color = :black, linestyle = :solid, subplot = 6, ylims = (0, maximum(decomp_pred.h_C)))
plot!(size = (600,300))
savefig("figures/"*folder*"siler_"*model*"_LEgrad_cC.pdf")


Lstar_plt = plot(layout = (1,2), legend = false,  margin=6Plots.mm)
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_c[pre_2020], yguidefontrotation=-90,
    xlabel = "Year", ylabel  = L"L^*_c", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_c[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (minimum(decomp_pred.Lstar_c), 0))
plot!(decomp_pred.year[pre_2020], decomp_pred.Lstar_C[pre_2020], yguidefontrotation=-90,
    xlabel = "Year", ylabel = L"L^*_C", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.Lstar_C[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0, maximum(decomp_pred.Lstar_C)))
plot!(size = (600,200))
savefig("figures/"*folder*"siler_"*model*"_Lstargrad.pdf")



# Plot the gradient of remaining LE at each age over time
all_years = decomp_pred.year
ages = 0:140
LEgrad_df = DataFrame(age = repeat(ages, length(all_years)), year = repeat(all_years, inner = length(ages)))
grad_vars = [:LE_Bs, :LE_bs, :LE_Cs, :LE_cs, :LE_ds,
    :Lstar_Bs, :Lstar_bs, :Lstar_Cs, :Lstar_cs, :Lstar_ds,
    :Lmed_Bs, :Lmed_bs, :Lmed_Cs, :Lmed_cs, :Lmed_ds,
    :h_Bs, :h_bs, :h_Cs, :h_cs, :h_ds,
    :H_Bs, :H_bs, :H_Cs, :H_cs, :H_ds]
LEgrad_df = hcat(LEgrad_df,DataFrame(NaN.*zeros(nrow(LEgrad_df), length(grad_vars)), grad_vars))
prog = Progress(length(all_years), desc = "Computing gradient df: ")
for yy in all_years
    # Get parameters for this year
    obs = Int.(1:nrow(decomp_pred))[decomp_pred.year .== yy][1]
    row = NamedTuple(decomp_pred[obs,:])
    params = SilerParam(b = row.b, B = row.B, c = row.c, C = row.C, d = row.d)

    # Fill grad dataframe
    out_obs = Int.(1:nrow(LEgrad_df))[LEgrad_df.year .== yy]
    # LE
    LEgrad_df.LE_Bs[out_obs] = LEgrad.([params], ages, [:B], spec = :Bergeron)
    LEgrad_df.LE_bs[out_obs] = LEgrad.([params], ages, [:b], spec = :Bergeron)
    LEgrad_df.LE_Cs[out_obs] = LEgrad.([params], ages, [:C], spec = :Bergeron)
    LEgrad_df.LE_cs[out_obs] = LEgrad.([params], ages, [:c], spec = :Bergeron)
    LEgrad_df.LE_ds[out_obs] = LEgrad.([params], ages, [:d], spec = :Bergeron)
    # LE
    LEgrad_df.Lstar_Bs[out_obs] .= lifespangrad.([params], [0.001], [:B], spec = :Bergeron)
    LEgrad_df.Lstar_bs[out_obs] .= lifespangrad.([params], [0.001], [:b], spec = :Bergeron)
    LEgrad_df.Lstar_Cs[out_obs] .= lifespangrad.([params], [0.001], [:C], spec = :Bergeron)
    LEgrad_df.Lstar_cs[out_obs] .= lifespangrad.([params], [0.001], [:c], spec = :Bergeron)
    LEgrad_df.Lstar_ds[out_obs] .= lifespangrad.([params], [0.001], [:d], spec = :Bergeron)
    # Lmed
    LEgrad_df.Lmed_Bs[out_obs] = Lmedgrad.([params], ages, [:B], spec = :Bergeron)
    LEgrad_df.Lmed_bs[out_obs] = Lmedgrad.([params], ages, [:b], spec = :Bergeron)
    LEgrad_df.Lmed_Cs[out_obs] = Lmedgrad.([params], ages, [:C], spec = :Bergeron)
    LEgrad_df.Lmed_cs[out_obs] = Lmedgrad.([params], ages, [:c], spec = :Bergeron)
    LEgrad_df.Lmed_ds[out_obs] = Lmedgrad.([params], ages, [:d], spec = :Bergeron)
    # H
    LEgrad_df.H_Bs[out_obs] = Hgrad.([params], ages, [:B], spec = :Bergeron)
    LEgrad_df.H_bs[out_obs] = Hgrad.([params], ages, [:b], spec = :Bergeron)
    LEgrad_df.H_Cs[out_obs] = Hgrad.([params], ages, [:C], spec = :Bergeron)
    LEgrad_df.H_cs[out_obs] = Hgrad.([params], ages, [:c], spec = :Bergeron)
    LEgrad_df.H_ds[out_obs] = Hgrad.([params], ages, [:d], spec = :Bergeron)
    # H
    LEgrad_df.h_Bs[out_obs] = hgrad.([params], ages, [:B], spec = :Bergeron)
    LEgrad_df.h_bs[out_obs] = hgrad.([params], ages, [:b], spec = :Bergeron)
    LEgrad_df.h_Cs[out_obs] = hgrad.([params], ages, [:C], spec = :Bergeron)
    LEgrad_df.h_cs[out_obs] = hgrad.([params], ages, [:c], spec = :Bergeron)
    LEgrad_df.h_ds[out_obs] = hgrad.([params], ages, [:d], spec = :Bergeron)
    next!(prog)
end

# C
plot(LEgrad_df.age, LEgrad_df.LE_Cs, group=LEgrad_df.year, legend = false,
    ylabel = L"LE_C(a)", xlabel = "Age")
hline!([0,0], color = :black, linestyle = :dash, label = false)

plot(LEgrad_df.age, LEgrad_df.LE_cs, group=LEgrad_df.year, legend = false,
    ylabel = L"LE_c(a)", xlabel = "Age")
hline!([0,0], color = :black, linestyle = :dash, label = false)

plot(LEgrad_df.age[LEgrad_df.age .<111], LEgrad_df.h_Cs[LEgrad_df.age .<111],
    group=LEgrad_df.year[LEgrad_df.age .<111], legend = false,
    ylabel = L"h_C(a)", xlabel = "Age")
hline!([0,0], color = :black, linestyle = :dash, label = false)

plot(LEgrad_df.age, LEgrad_df.h_cs, group=LEgrad_df.year, legend = false,
    ylabel = L"h_c(a)", xlabel = "Age")
hline!([0,0], color = :black, linestyle = :dash, label = false)


CSV.write("figures/"*folder*"/siler_"*model*"_LEgrads.csv", LEgrad_df)




"""
Cross derivatives
"""
# Initialise the parameters at the latest estimated level
obs = decomp_pred.year .== 2018
param_base = SilerParam(b = decomp_pred.b[obs][1], B = decomp_pred.B[obs][1],
    c = decomp_pred.c[obs][1], C = decomp_pred.C[obs][1], d = decomp_pred.d[obs][1])


param = deepcopy(param_base)
plot(ages, LEcross.([param], ages, [:c], [:C], spec = :Bergeron), xlabel = "Age",
 ylabel = L"LE_{cC}(a)")

plot(ages, repeat([0.], length(ages)), xlabel = "Age", linestyle = :dash, color = :black,
    ylabel = L"LE_{C}(a)", label = false)
for c_ in 0.07:0.01:0.14
    param = deepcopy(param_base)
    param.c = c_
    plot!(ages, LEgrad.([param], ages, [:C], spec = :Bergeron), label = "c = "*string(param.c))
end
plot!(title = "Gradient wrt C for different c")

plot(ages, repeat([0.], length(ages)), xlabel = "Age", linestyle = :dash, color = :black,
    ylabel = L"LE_{c}(a)", label = false)
for C_ in 80:3.0:110.0
    param = deepcopy(param_base)
    param.C = C_
    plot!(ages, LEgrad.([param], ages, [:c], spec = :Bergeron), label = "C = "*string(param.C))
end
plot!(title = "Gradient wrt c for different C")



plot(ages, siler_S.([param], [0.0], ages, spec = :Bergeron))




"""
Out-of-sample forecasts
"""
oos_forecasts = false
if oos_forecasts
    folder = "benchmark/held-out"
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
        niters = 1250
        nchains = 4
        # MCMC sampling
        chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
            niters, nchains, init_params = prior_i2_vals)
        display(chain_i2)
        @save "figures/"*folder*"/siler_"*model*"_chain_"*string(years_selected[T-held_out])*".jld2" chain_i2

    end
end


"""
End of script
"""
