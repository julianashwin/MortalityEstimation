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


"""
Import data
"""
# Which data do we want to look at
code = "Best Practice"
folder = "figures/benchmark/"
model = "indep"
# Full raw samples
@load folder*"siler_"*model*"_chain.jld2" chain_indep

# Some more convenient summary statistics for the sampled posterior
parests_df_ber =  CSV.read(folder*"siler_"*model*"_params_ber.csv", DataFrame, ntasks = 1)
parests_df_col =  CSV.read(folder*"siler_"*model*"_params_col.csv", DataFrame, ntasks = 1)
parests_df_sco =  CSV.read(folder*"siler_"*model*"_params_sco.csv", DataFrame, ntasks = 1)
# Also the underlying mortality data for comparison
mort_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1)
mort_df = mort_df[mort_df.best_practice.==1,:]
sort!(mort_df, [:year, :age])

"""
Create decomposition dataframes for each specification
"""
## Parameters to look at infant mortality
decomp_df_col = create_decomp(parests_df_col; spec = :Colchero, eval_age = 0)
decomp_df_sco = create_decomp(parests_df_sco; spec = :Scott, eval_age = 0)
decomp_df_ber = create_decomp(parests_df_ber; spec = :Bergeron, eval_age = 0)
# Plot model implied LE as a sense check (should be the same across as specifications
scatter(decomp_df_col.year, decomp_df_col.LE_mod, markershape = :circle)
scatter!(decomp_df_sco.year, decomp_df_sco.LE_mod, markershape = :cross)
scatter!(decomp_df_ber.year, decomp_df_ber.LE_mod, markershape = :xcross)

"""
Plot the historical decompositions
"""
# Bergeron
le_p = plot_decomp(decomp_df_ber, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_ber, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Bergeron parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig(folder*"siler_"*model*"_decomp_ber.pdf")
# Colchero
le_p = plot_decomp(decomp_df_col, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_col, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Colchero parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig(folder*"siler_"*model*"_decomp_col.pdf")
# Scott
le_p = plot_decomp(decomp_df_sco, :LE)
le_p = plot!(legend = false, title = "Life Expectancy")
h_p = plot_decomp(decomp_df_sco, :H)
h_p = plot!(title = "Lifespan Inequality")
p_title = plot(title = "Historical decomposition Siler Scott parameters "*string(code),
    grid = false, showaxis = false, bottom_margin = -10Plots.px, yticks = false, xticks = false)
display(plot(p_title, le_p, h_p, layout = @layout([A{0.01h}; B C]), size = (800,400)))
savefig(folder*"siler_"*model*"_decomp_sco.pdf")


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
# Bergeron
le_p = plot_LEgrads(decomp_df_ber)
savefig(folder*"siler_"*model*"_LEgrad_ber.pdf")
# Colchero
le_p = plot_LEgrads(decomp_df_ber)
savefig(folder*"siler_"*model*"_LEgrad_col.pdf")
# Scott
le_p = plot_LEgrads(decomp_df_sco)
savefig(folder*"siler_"*model*"_LEgrad_sco.pdf")
