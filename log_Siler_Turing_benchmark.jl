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
#chain_i2_cov = deepcopy(chain_i2)
#@load "figures/"*folder*"/siler_i2_chain.jld2" chain_i2
#chain_i2
folder = "benchmark"

plot(chain_i2[["σ_pars[3]","σ_pars[4]", "ρ_cC", "ρ_cd", "ρ_Cd"]])
df_post_cov = DataFrame(chain_i2)
#df_post_cov = DataFrame(chain_i2_cov)
#cor(df_post[df_post.chain .== 3,"lB[5]"],df_post[df_post.chain.==3,"lb[5]"])
#cor(df_post[df_post.chain .== 4,"lc[5]"],df_post[df_post.chain.==4,"ld[5]"])


cor(df_post_cov[df_post_cov.chain .== 4,"lpars[5][1]"],df_post_cov[df_post_cov.chain.==4,"lpars[5][2]"])
cor(df_post_cov[df_post_cov.chain .== 1,"lpars[4][3]"],df_post_cov[df_post_cov.chain.==1,"lpars[4][4]"])
mean(df_post_cov[df_post_cov.chain .>= 0,"ρ_cC"])
cor(df_post_cov[df_post_cov.chain .== 1,"α_C[4]"],df_post_cov[df_post_cov.chain.==1,"α_c[4]"])

"""
Plot Siler parameters
"""
# Display the parameters with Bergeron specification
parests_i2_ber = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :cov, spec = :Bergeron)
p1 = plot_siler_params(parests_i2_ber)
p_title = plot(title = "I(2) Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_ber.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_ber.csv", parests_i2_ber)
# with Colchero specification
parests_i2_col = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :cov, spec = :Colchero)
p1 = plot_siler_params(parests_i2_col)
p_title = plot(title = "I(2) Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_"*model*"_params_col.pdf")
CSV.write("figures/"*folder*"/siler_"*model*"_params_col.csv", parests_i2_col)


## Plot time series parameters
p2 = plot_ts_params(parests_i2_ber, model_vers = :cov)
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
df_post = compute_LE_post(df_post, years_selected, nahead, spec = :Bergeron, model_vers = :cov)
# Compute df_pred which extends df_post with predictions
df_pred = compute_forecasts(df_post, nahead, ndraws, years_selected, spec = :Bergeron, model_vers = :cov)
# Summarize past and forecast variables
past_years = years_selected
fut_years = Int.(maximum(years_selected) .+ 5.0.*(1:nahead))
parests_pred = extract_forecast_variables(df_pred, past_years, fut_years,
    log_pars = true, spec = :Bergeron, model_vers = :cov)
CSV.write("figures/"*folder*"/siler_"*model*"_preds.csv", parests_pred)

# Plot Siler parameter forecasts
plot_siler_params(parests_pred, forecasts = true)
savefig("figures/"*folder*"/siler_"*model*"_param_pred.pdf")
# Plot time series parameter forecasts
plot_ts_params(parests_pred, model_vers = :cov, forecasts = true)
savefig("figures/"*folder*"/siler_"*model*"_ts_pred.pdf")
# Plot forecasts for model implied LE and H
plot_LE_H(parests_pred, forecasts = true, bands = true)
savefig("figures/"*folder*"/siler_"*model*"_leh_pred.pdf")
# Forecast decomposition of future LE and H
decomp_pred = create_decomp(parests_pred[parests_pred.year .> 1985,:],
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



"""
LE gradients over time
"""
# prediction decomposition
decomp_pred = create_decomp(parests_pred, spec = :Bergeron, eval_age = 0, forecasts = true)
pre_2020 = decomp_pred.year .<= 2018
post_2020 = decomp_pred.year .>= 2018
# Bergeron
bB_plt = plot(layout = (2,2), legend = false)
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_b[pre_2020], title = L"LE_b", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_b[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,maximum(decomp_pred.LE_b)))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_B[pre_2020], title = L"LE_B", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_B[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,maximum(decomp_pred.LE_B)))
plot!(decomp_pred.year[pre_2020], decomp_pred.H_b[pre_2020], title = L"H_b", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.H_b[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, ylims = (minimum(decomp_pred.H_b),0), xmirror = true)
plot!(decomp_pred.year[pre_2020], decomp_pred.H_B[pre_2020], title = L"H_B", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.H_B[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, ylims = (minimum(decomp_pred.H_B),0), xmirror = true)
plot!(size = (600,300))
savefig("figures/"*folder*"siler_"*model*"_LEgrad_bB.pdf")


cC_plt = plot(layout = (2,2), legend = false)
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_c[pre_2020], title = L"LE_c", subplot = 1)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_c[post_2020], linestyle = :dash,subplot = 1)
hline!([0,0], color = :black, linestyle = :solid, subplot = 1, ylims = (0,75))
plot!(decomp_pred.year[pre_2020], decomp_pred.LE_C[pre_2020], title = L"LE_C", subplot = 2)
plot!(decomp_pred.year[post_2020], decomp_pred.LE_C[post_2020], linestyle = :dash,subplot = 2)
hline!([0,0], color = :black, linestyle = :solid, subplot = 2, ylims = (0,1))
plot!(decomp_pred.year[pre_2020], decomp_pred.h_c[pre_2020], title = L"h_c", subplot = 3)
plot!(decomp_pred.year[post_2020], decomp_pred.h_c[post_2020], linestyle = :dash,subplot = 3)
hline!([0,0], color = :black, linestyle = :solid, subplot = 3, xmirror = true)#, ylims = (-2.4,0))
plot!(decomp_pred.year[pre_2020], decomp_pred.h_C[pre_2020], title = L"h_C", subplot = 4)
plot!(decomp_pred.year[post_2020], decomp_pred.h_C[post_2020], linestyle = :dash,subplot = 4)
hline!([0,0], color = :black, linestyle = :solid, subplot = 4, xmirror = true)#, ylims = (-0.002,0),
plot!(size = (600,300))
savefig("figures/"*folder*"siler_"*model*"_LEgrad_cC.pdf")


# Plot the gradient of remaining LE at each age over time
all_years = decomp_pred.year
ages = 0:110
LEgrad_df = DataFrame(age = repeat(ages, length(all_years)), year = repeat(all_years, inner = length(ages)))
grad_vars = [:LE_Bs, :LE_bs, :LE_Cs, :LE_cs, :LE_ds]
LEgrad_df = hcat(LEgrad_df,DataFrame(NaN.*zeros(nrow(LEgrad_df), length(grad_vars)), grad_vars))
for yy in all_years
    # Get parameters for this year
    obs = Int.(1:nrow(decomp_pred))[decomp_pred.year .== yy][1]
    row = NamedTuple(decomp_pred[obs,:])
    params = SilerParam(b = row.b, B = row.B, c = row.c, C = row.C, d = row.d)

    # Fill grad dataframe
    out_obs = Int.(1:nrow(LEgrad_df))[LEgrad_df.year .== yy]
    LEgrad_df.LE_Bs[out_obs] = LEgrad.([params], ages, [:B], spec = :Bergeron)
    LEgrad_df.LE_bs[out_obs] = LEgrad.([params], ages, [:b], spec = :Bergeron)
    LEgrad_df.LE_Cs[out_obs] = LEgrad.([params], ages, [:C], spec = :Bergeron)
    LEgrad_df.LE_cs[out_obs] = LEgrad.([params], ages, [:c], spec = :Bergeron)
    LEgrad_df.LE_ds[out_obs] = LEgrad.([params], ages, [:d], spec = :Bergeron)
end

plot(LEgrad_df.age, LEgrad_df.LE_Cs, group=LEgrad_df.year)
plot(LEgrad_df.age, LEgrad_df.LE_cs, group=LEgrad_df.year)

CSV.write("figures/"*folder*"/siler_"*model*"_LEgrads.csv", LEgrad_df)


"""
End of script
"""
