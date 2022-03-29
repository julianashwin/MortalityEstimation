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
Plot priors
"""
## Set priors
# Mean of IG is b/(a-1)
plot(layout = (2,3), yticks = false, size = (1000,500))
plot!(LogNormal(log(5), 2.0), title = L"B_{1} \sim LogNormal(ln(5),2)",
    label = false, subplot = 1, xlims = (0,100))
plot!(LogNormal(log(1), 1.0), xlim=(0,10), title = L"b_{1} \sim LogNormal(ln(1),2)",
    label = false, subplot = 2)
plot!(LogNormal(log(10), 2.0), title = L"C_{1} \sim LogNormal(ln(120),2)",
    label = false, subplot = 4, xlims = (0,200))
plot!(LogNormal(log(0.1), 1.0), xlim=(0,1), title = L"c_{1} \sim LogNormal(ln(0.1),2)",
    label = false, subplot = 5)
plot!(LogNormal(log(0.025), 1.0), title = L"d_{1} \sim LogNormal(ln(0.025),2)",
    label = false, subplot = 3)
plot!(LogNormal(log(0.001), 1.0), title = L"\sigma_{1} \sim LogNormal(ln(0.001),2)",
    label = false, subplot = 6)
savefig("figures/general/log_siler_priors.pdf")


"""
Set up a single country/case to run through examples
"""
## Data prep for single coujntry
code = "Best Practice"
folder = "best_practice"
#country_df = mort_df[(mort_df.code .== code), :]
country_df = bp_alt_df
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
Static Siler model
"""
# Choose year
yy = 1
## Find some starting points
# MAP estimate for static model on raw mortality
map_static = optimize(siler_static((country_m_data[yy]), country_ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
map_vals =  map_static.values.array
map_param = SilerParam(b = map_vals[2], B = map_vals[1], c = map_vals[4], C = map_vals[3],
    d = map_vals[5], σ = map_vals[6])
# Map estimate for static model on log mortality
map_log_static = optimize(log_siler_static((country_lm_data[yy]), country_ages), MAP(),
    Optim.Options(iterations=50_000, allow_f_increases=true))
map_log_vals =  map_log_static.values.array
map_log_param = SilerParam(b = map_log_vals[2], B = map_log_vals[1], c = map_log_vals[4], C = map_log_vals[3],
    d = map_log_vals[5], σ = map_log_vals[6])
# Alternatively, we can simulate from the prior and start there
prior_static = sample(siler_static(country_m_data[yy], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_static)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_vals = df_prior[1,3:(end-1)]
prior_param = SilerParam(b = prior_vals.b, B = prior_vals.B, c = prior_vals.c, C = prior_vals.C,
    d = prior_vals.d, σ = prior_vals.σ)


## Estimate models
# Number of MCMC iterations and chains
niters = 1000
nchains = 1
# Raw mortality model
chain_static = sample(siler_static(country_m_data[1], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = map_log_vals)
df_chain = DataFrame(chain_static)
insert!.(eachcol(df_chain), 1, vcat([0,0],median.(eachcol(df_chain[:,3:end]))))
chain_vals = df_chain[1,:]
chain_param = SilerParam(b = chain_vals.b, B = chain_vals.B, c = chain_vals.c, C = chain_vals.C,
    d = chain_vals.d, σ = chain_vals.σ)
# Raw mortality model
chain_log_static = sample(log_siler_static(country_lm_data[1], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = map_log_vals)
df_chain = DataFrame(chain_log_static)
insert!.(eachcol(df_chain), 1, vcat([0,0],median.(eachcol(df_chain[:,3:end]))))
chain_log_vals = df_chain[1,:]
chain_log_param = SilerParam(b = chain_log_vals.b, B = chain_log_vals.B, c = chain_log_vals.c, C = chain_log_vals.C,
    d = chain_log_vals.d, σ = chain_log_vals.σ)

## Plot fit
scatter(country_m_data[yy], markershape = :cross, markeralpha = 0.5,
    label = "Data ("*string(country_years[1])*")", legend = :topleft)
plot!(siler.([prior_param], 0:110, spec = :Colchero), label = "Static model prior")
plot!(siler.([map_param], 0:110, spec = :Colchero), label = "Static model MAP", linestyle = :dash)
plot!(siler.([chain_param], 0:110, spec = :Colchero), label = "Static model post. median", linestyle = :dot)
plot!(siler.([map_log_param], 0:110, spec = :Colchero), label = "Static log model MAP", linestyle = :dash)
plot!(siler.([chain_log_param], 0:110, spec = :Colchero), label = "Static log model post. median", linestyle = :dot)
savefig("figures/"*folder*"/siler_vs_log_siler_fit.pdf")


"""
Multiple independent Siler models
"""
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
    niters, nchains, init_params = map_indep_vals)
display(chain_indep)
CSV.write("figures/"*folder*"/siler_indep_fullpost.csv", DataFrame(chain_indep))
## Plot Siler parameters
# Display the parameters with Colchero specification
parests_indep_col = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Colchero)
p1 = plot_siler_params(parests_indep_col)
p_title = plot(title = "Multiple independent Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_indep_params_col.pdf")
CSV.write("figures/"*folder*"/siler_indep_params_col.csv", parests_indep_col)
# With Scott specification
parests_indep_sco = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Scott)
p1 = plot_siler_params(parests_indep_sco)
p_title = plot(title = "Multiple independent Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_indep_params_sco.pdf")
CSV.write("figures/"*folder*"/siler_indep_params_sco.csv", parests_indep_sco)
# With Bergeron specification
parests_indep_ber = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Bergeron)
p1 = plot_siler_params(parests_indep_ber)
p_title = plot(title = "Multiple independent Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_indep_params_ber.pdf")
CSV.write("figures/"*folder*"/siler_indep_params_ber.csv", parests_indep_ber)
# With Standard specification
parests_indep_sta = extract_variables(chain_indep, years_selected, log_pars = false,
    model_vers = :indep, spec = :Standard)
p1 = plot_siler_params(parests_indep_sta)
p_title = plot(title = "Multiple independent Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_indep_params_sta.pdf")
CSV.write("figures/"*folder*"/siler_indep_params_sta.csv", parests_indep_sta)




"""
Dynamic Siler model with simple random walk
"""
# Adjust if you don't want every period
periods = Int.(1:7:T)
years_selected = Int.(round.(country_years[periods]))
## Find some starting points
# MAP estimate for multiple independent models on log mortality
map_rw = optimize(log_siler_justrw(country_lm_data[periods], country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=60_000, allow_f_increases=true))
map_rw_vals =  map_rw.values.array
# Alternatively, we can simulate from the prior and start there
prior_rw = sample(log_siler_justrw(country_lm_data[periods], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_rw)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_rw_vals = df_prior[1,3:(end-1)]
## Estimate the model
niters = 1000
nchains = 4
# MCMC sampling
chain_rw = sample(log_siler_justrw(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_rw_vals)
display(chain_rw)

## Plot Siler parameters
# Display the parameters with Colchero specification
parests_rw_col = extract_variables(chain_rw, years_selected, log_pars = true,
    model_vers = :justrw, spec = :Colchero)
p1 = plot_siler_params(parests_indep_col)
p_title = plot(title = "Random walk Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_rw_params_col.pdf")
# With Scott specification
parests_rw_sco = extract_variables(chain_rw, years_selected, log_pars = true,
    model_vers = :justrw, spec = :Scott)
p1 = plot_siler_params(parests_rw_sco)
p_title = plot(title = "Random walk Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_rw_params_sco.pdf")
# With Bergeron specification
parests_rw_ber = extract_variables(chain_rw, years_selected, log_pars = true,
    model_vers = :justrw, spec = :Bergeron)
p1 = plot_siler_params(parests_rw_ber)
p_title = plot(title = "Random walk Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_rw_params_ber.pdf")
# With Standard specification
parests_rw_sta = extract_variables(chain_rw, years_selected, log_pars = true,
    model_vers = :justrw, spec = :Standard)
p1 = plot_siler_params(parests_rw_sta)
p_title = plot(title = "Random walk Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_rw_params_sta.pdf")

## Plot time series parameters
p2 = plot_ts_params(parests_rw_col, model_vers = :justrw)
p_title = plot(title = "Random walk Siler Colchero ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_rw_ts_params_col.pdf")









"""
Dynamic Siler model with parameters in first differences
"""
# Adjust if you don't want every period
periods = Int.(1:7:T)
years_selected = Int.(round.(country_years[periods]))
## Find some starting points
# MAP estimate for multiple independent models on log mortality
map_fd = optimize(log_siler_dyn_firstdiff(country_lm_data[periods], country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=60_000, allow_f_increases=true))
map_fd_vals =  map_fd.values.array
# Alternatively, we can simulate from the prior and start there
prior_fd = sample(log_siler_dyn_firstdiff(country_lm_data[periods], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_fd)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_fd_vals = df_prior[1,3:(end-1)]
## Estimate the model
niters = 1000
nchains = 4
# MCMC sampling
chain_fd = sample(log_siler_dyn_firstdiff(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_fd_vals)
display(chain_fd)

## Plot Siler parameters
# Display the parameters with Colchero specification
parests_fd_col = extract_variables(chain_fd, years_selected, log_pars = true,
    model_vers = :firstdiff, spec = :Colchero)
p1 = plot_siler_params(parests_fd_col)
p_title = plot(title = "First differences Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_fd_params_col.pdf")
# With Scott specification
parests_fd_sco = extract_variables(chain_fd, years_selected, log_pars = true,
    model_vers = :firstdiff, spec = :Scott)
p1 = plot_siler_params(parests_fd_sco)
p_title = plot(title = "First differences Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_fd_params_sco.pdf")
# With Bergeron specification
parests_fd_ber = extract_variables(chain_fd, years_selected, log_pars = true,
    model_vers = :firstdiff, spec = :Bergeron)
p1 = plot_siler_params(parests_fd_ber)
p_title = plot(title = "First differences Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_fd_params_ber.pdf")
# With Standard specification
parests_fd_sta = extract_variables(chain_fd, years_selected, log_pars = true,
    model_vers = :firstdiff, spec = :Standard)
p1 = plot_siler_params(parests_fd_sta)
p_title = plot(title = "First differences Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_fd_params_sta.pdf")

## Plot time series parameters
p2 = plot_ts_params(parests_fd_col, model_vers = :firstdiff)
p_title = plot(title = "First differences Siler Colchero ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_fd_ts_params_col.pdf")



















"""
Dynamic Siler model with parameters i(2) random walk with drift
"""
# Adjust if you don't want every period
periods = Int.(1:T)
years_selected = Int.(round.(country_years[periods]))
## Find some starting points
# MAP estimate for multiple independent models on log mortality
map_i2 = optimize(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), MAP(), LBFGS(),
    Optim.Options(iterations=60_000, allow_f_increases=true))
map_i2_vals =  map_i2.values.array
# Alternatively, we can simulate from the prior and start there
prior_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), Prior(), 5000)
df_prior = DataFrame(prior_i2)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_i2_vals = df_prior[1,3:(end-1)]
## Estimate the model
niters = 1000
nchains = 4
# MCMC sampling
chain_i2 = sample(log_siler_dyn_i2drift(country_lm_data[periods], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains, init_params = prior_i2_vals)
display(chain_i2)
CSV.write("figures/"*folder*"/siler_i2_fullpost.csv", DataFrame(chain_i2))

## Plot Siler parameters
# Display the parameters with Colchero specification
parests_i2_col = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
p1 = plot_siler_params(parests_i2_col)
p_title = plot(title = "I(2) Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_i2_params_col.pdf")
CSV.write("figures/"*folder*"/siler_i2_params_col.csv", parests_indep_sco)
# With Scott specification
parests_i2_sco = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Scott)
p1 = plot_siler_params(parests_i2_sco)
p_title = plot(title = "I(2) Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_i2_params_sco.pdf")
CSV.write("figures/"*folder*"/siler_i2_params_sco.csv", parests_indep_sco)
# With Bergeron specification
parests_i2_ber = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Bergeron)
p1 = plot_siler_params(parests_i2_ber)
p_title = plot(title = "I(2) Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_i2_params_ber.pdf")
CSV.write("figures/"*folder*"/siler_i2_params_ber.csv", parests_indep_sco)
# With Standard specification
parests_i2_sta = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Standard)
p1 = plot_siler_params(parests_i2_sta)
p_title = plot(title = "I(2) Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_i2_params_sta.pdf")
CSV.write("figures/"*folder*"/siler_i2_params_sta.csv", parests_indep_sco)

## Plot time series parameters
p2 = plot_ts_params(parests_i2_col, model_vers = :i2drift)
p_title = plot(title = "I(2) Siler Colcher ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/"*folder*"/siler_i2_ts_params_col.pdf")








"""
End of script
"""
