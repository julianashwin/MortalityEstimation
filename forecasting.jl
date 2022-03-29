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
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings, JLD2

include("src/MortalityEstimation.jl")
gr()

code = "Best Practice"
folder = "figures/best_practice/"
model = "i2"

# Import the full posterior
@load folder*"siler_"*model*"_chain.jld2" chain_i2
chain_post = deepcopy(chain_i2)
parests_df_col =  CSV.read(folder*"siler_"*model*"_params_col.csv", DataFrame, ntasks = 1)
#parests_df_col.year[parests_df_col.year .== 0] .= missing

years = unique(parests_df_col.year[parests_df_col.year .> 0])

"""
Identify the variable names
"""
## Siler parameters
Bs = Symbol.("lB[".*string.(1:length(years)).*"]")
bs = Symbol.("lb[".*string.(1:length(years)).*"]")
Cs = Symbol.("lC[".*string.(1:length(years)).*"]")
cs = Symbol.("lc[".*string.(1:length(years)).*"]")
ds = Symbol.("ld[".*string.(1:length(years)).*"]")
σs = Symbol.("lσ[".*string.(1:length(years)).*"]")
## Drift terms
α_Bs = Symbol.("α_B[".*string.(1:length(years)-1).*"]")
α_bs = Symbol.("α_b[".*string.(1:length(years)-1).*"]")
α_Cs = Symbol.("α_C[".*string.(1:length(years)-1).*"]")
α_cs = Symbol.("α_c[".*string.(1:length(years)-1).*"]")
α_ds = Symbol.("α_d[".*string.(1:length(years)-1).*"]")
α_σs = Symbol.("α_σ[".*string.(1:length(years)-1).*"]")
## Variance terms
σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
σ_α_par_names = Symbol.("σ_α".*["B", "b", "C", "c", "d", "σ"])


"""

"""



σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
    Int.(1:6), parname = :σ_par)
σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
σ_pars_ests.year .= 0

par_ests = vcat(par_ests, σ_pars_ests)

σ_α_pars_ests = summarise_stats(df_post[:,(σ_α_par_names)], stats[in.(stats.parameters, [σ_α_par_names]),:],
    Int.(1:6), parname = :σ_α_par)
σ_α_pars_ests.parameter = [:σ_αB, :σ_αb, :σ_αC, :σ_αc, :σ_αd, :σ_ασ]
σ_α_pars_ests.year .= 0

par_ests = vcat(par_ests, σ_α_pars_ests)


α_B_ests = summarise_stats(df_post[:,α_Bs], stats[in.(stats.parameters, [α_Bs]),:],
    years[2:end], parname = :α_B)
α_b_ests = summarise_stats(df_post[:,α_bs], stats[in.(stats.parameters, [α_bs]),:],
    years[2:end], parname = :α_b)
α_C_ests = summarise_stats(df_post[:,α_Cs], stats[in.(stats.parameters, [α_Cs]),:],
    years[2:end], parname = :α_C)
α_c_ests = summarise_stats(df_post[:,α_cs], stats[in.(stats.parameters, [α_cs]),:],
    years[2:end], parname = :α_c)
α_d_ests = summarise_stats(df_post[:,α_ds], stats[in.(stats.parameters, [α_ds]),:],
    years[2:end], parname = :α_d)
α_σ_ests = summarise_stats(df_post[:,α_σs], stats[in.(stats.parameters, [α_σs]),:],
    years[2:end], parname = :α_σ)

par_ests = vcat(par_ests, α_B_ests, α_b_ests, α_C_ests, α_c_ests, α_d_ests, α_σ_ests)

samples_df[:, Bs] = exp.(samples_df[:, Bs])
samples_df[:, bs] = exp.(samples_df[:, bs])
samples_df[:, Cs] = exp.(samples_df[:, Cs])
samples_df[:, cs] = exp.(samples_df[:, cs])
samples_df[:, ds] = exp.(samples_df[:, ds])
samples_df[:, σs] = exp.(samples_df[:, σs])
