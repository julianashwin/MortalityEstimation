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
all_df = CSV.read("data/clean/all_lifetab_1y.csv", DataFrame, ntasks = 1)

# Best practice data
mort_df = all_df[(all_df.year .>= 1900), :]
sort!(mort_df, [:code, :year, :age])
bp_df = mort_df[(mort_df.best_practice .== 1), :]

sort!(bp_df, [:year, :age])

"""
Toggle options
"""
# Choose here what data to use
code = "Best Practice"
select_df = bp_df
m_var = :mx_f

plot(select_df.age, select_df[:,m_var], group = select_df.year, legend = :top)

"""
Preliminary checks and illustrative estimation of first and last period
"""
## Data prep for single country
# Suprisingly, we actually have some zeros here for small countries (e.g. ISL)
select_df[select_df[:,m_var] .== 0.0,m_var] .=  minimum(select_df[select_df[:,m_var] .> 0.0,m_var])
# Convert this into a matrix of mortality rates over time, age and year vectors
m_data = chunk(select_df[:,m_var], 110)
plot(m_data, legend = false)
lm_data = [log.(m_dist) for m_dist in m_data]
plot(lm_data, legend = false)
ages = Int64.(0:maximum(select_df.age))
years = unique(select_df.year)
T = length(lm_data)

@assert length(m_data)==length(years) "number of years doesn't match length of m_data"
@assert length(m_data[1])==length(ages) "number of ages doesn't match length of m_data[1]"

## Estimate model
@time map_est = optimize(log_siler_dyn_i2drift(lm_data, ages), MAP(), LBFGS(),
    Optim.Options(iterations=60_000, allow_f_increases=true))
print("Estimated MAP for "*code)
map_vals =  map_est.values.array
showtable(DataFrame(name = names(coef(map_est))[1], value = coef(map_est)))

prior_est = sample(log_siler_dyn_i2drift(lm_data, ages), Prior(), 5000)
df_prior = DataFrame(prior_est)
insert!.(eachcol(df_prior), 1, vcat([0,0],median.(eachcol(df_prior[:,3:end]))))
prior_vals = df_prior[1,3:(end-1)]
# Estimate by MCMC
@time chain_i2 = sample(log_siler_dyn_i2drift(lm_data, ages), NUTS(0.65), MCMCThreads(),
    1000, 4, init_params = prior_vals)
display(chain_i2)
CSV.write("figures/best_practice/siler_i2_params_fullpost.csv", DataFrame(chain_i2))
years_selected = country_years

## Plots the Siler parameters
# Display the parameters with Colchero specification
parests_i2_col = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Colchero)
p1 = plot_siler_params(parests_i2_col)
p_title = plot(title = "I(2) Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/best_practice/siler_i2_params_col.pdf")
CSV.write("figures/best_practice/siler_i2_params_col.csv", parests_i2_col)
# With Scott specification
parests_i2_sco = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Scott)
p1 = plot_siler_params(parests_i2_sco)
p_title = plot(title = "I(2) Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/best_practice/siler_i2_params_sco.pdf")
CSV.write("figures/best_practice/siler_i2_params_sco.csv", parests_i2_sco)
# With Bergeron specification
parests_i2_ber = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Bergeron)
p1 = plot_siler_params(parests_i2_ber)
p_title = plot(title = "I(2) Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/best_practice/siler_i2_params_ber.pdf")
CSV.write("figures/best_practice/siler_i2_params_ber.csv", parests_i2_ber)
# With Standard specification
parests_i2_sta = extract_variables(chain_i2, years_selected, log_pars = true,
    model_vers = :i2drift, spec = :Standard)
p1 = plot_siler_params(parests_i2_sta)
p_title = plot(title = "I(2) Siler Standard parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p1, layout = @layout([A{0.01h}; B])))
savefig("figures/best_practice/siler_i2_params_sta.pdf")
CSV.write("figures/best_practice/siler_i2_params_sta.csv", parests_i2_sta)

## Plot time series parameters
p2 = plot_ts_params(parests_i2_col, model_vers = :i2drift)
p_title = plot(title = "I(2) Siler Colcher ts params "*string(code), grid = false, showaxis = false,
    bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, p2, layout = @layout([A{0.01h}; B])))
savefig("figures/best_practice/siler_i2_ts_params_col.pdf")


## Plot a decomposition of LE and H



## Compute some forecasts





"""
End of script
"""
