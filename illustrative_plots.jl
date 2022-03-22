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
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings

include("src/MortalityEstimation.jl")
gr()

"""
Select parameters
"""
## Parameters to look at infant mortality
param_col_1 = SilerParam(b = 1.0, B = 1.25, c = 0.15, C = 10.0, d = 0.001)
param_sco_1 = SilerParam(b = param_col_1.b, B = param_col_1.B/param_col_1.b,
    c = param_col_1.c, C = param_col_1.C/param_col_1.c, d = param_col_1.d)
param_ber_1 = SilerParam(b = param_col_1.b, B= exp(-param_col_1.B),
    c = param_col_1.c, C = (param_col_1.C + log(param_col_1.c))/param_col_1.c,
    d = param_col_1.d)
plot(siler.([param_col_1], 0:110, spec = :Colchero), label = "Colchero")
plot!(siler.([param_sco_1], 0:110, spec = :Scott), label = "Scott")
plot!(siler.([param_ber_1], 0:110, spec = :Bergeron), label = "Bergeron")
LE_base_1 = LE(param_col_1, 0.0, spec = :Colchero)
H_base_1 = H(param_col_1, 0., spec = :Colchero)

## Parameters to look at senescent mortality
param_col_2 = SilerParam(b = 2.0, B = 3.5, c = 0.1, C = 10.0, d = 0.001)
param_sco_2 = SilerParam(b = param_col_2.b, B = param_col_2.B/param_col_2.b,
    c = param_col_2.c, C = param_col_2.C/param_col_2.c, d = param_col_2.d)
param_ber_2 = SilerParam(b = param_col_2.b, B= exp(-param_col_2.B),
    c = param_col_2.c, C = (param_col_2.C + log(param_col_2.c))/param_col_2.c,
    d = param_col_2.d)
plot(siler.([param_col_2], 0:110, spec = :Colchero), label = "Colchero")
plot!(siler.([param_sco_2], 0:110, spec = :Scott), label = "Scott")
plot!(siler.([param_ber_2], 0:110, spec = :Bergeron), label = "Bergeron")
LE_base_2 = LE(param_col_2, 0.0, spec = :Colchero)
H_base_2 = H(param_col_2, 0., spec = :Colchero)


"""
Initialise baseline for each case
"""
## Initialise dataframe for infant mortality
col_df_1 = init_illus(param_col_1, :Colchero)
base_plt(col_df_1[1:31,:], param_col_1); plot!(title = "Colchero infant baseline")
sco_df_1 = init_illus(param_sco_1, :Scott)
base_plt(sco_df_1[1:31,:], param_sco_1); plot!(title = "Scott infant baseline")
ber_df_1 = init_illus(param_ber_1, :Bergeron)
base_plt(ber_df_1[1:31,:], param_ber_1); plot!(title = "Bergeron infant baseline")

## Initialise dataframe for senescent mortality
col_df_2 = init_illus(param_col_2, :Colchero)
base_plt(col_df_2, param_col_2); plot!(title = "Colchero senescent baseline")
sco_df_2 = init_illus(param_sco_2, :Scott)
base_plt(sco_df_2, param_sco_2); plot!(title = "Scott senescent baseline")
ber_df_2 = init_illus(param_ber_2, :Bergeron)
base_plt(ber_df_2, param_ber_2); plot!(title = "Bergeron senescent baseline")



"""
Colchero specification example
"""
# Compute parameter changes
col_df_1, b_plt = param_change_plt(col_df_1[1:31,:], param_col_1, :b; LE_change = 5, spec = :Colchero)
col_df_1, B_plt = param_change_plt(col_df_1, param_col_1, :B; LE_change = 5, spec = :Colchero)
col_df_2, c_plt = param_change_plt(col_df_2, param_col_2, :c; LE_change = 10, spec = :Colchero)
col_df_2, C_plt = param_change_plt(col_df_2, param_col_2, :C; LE_change = 10, spec = :Colchero)
# Plot
l = @layout [a{0.48w} _ b{0.48w}; c{0.48w} _ d{0.48w}]
plot(b_plt, c_plt, B_plt, C_plt, layout  = l, size = (1200,800), margin = 40Plots.mm)
savefig("figures/interpret/colchero_example.pdf")


"""
Scott specification example
"""
# Compute parameter changes
sco_df_1, b_plt = param_change_plt(sco_df_1[1:31,:], param_sco_1, :b; LE_change = 5, spec = :Scott)
sco_df_1, B_plt = param_change_plt(sco_df_1, param_sco_1, :B; LE_change = 5, spec = :Scott)
sco_df_2, c_plt = param_change_plt(sco_df_2, param_sco_2, :c; LE_change = 10, spec = :Scott)
sco_df_2, C_plt = param_change_plt(sco_df_2, param_sco_2, :C; LE_change = 10, spec = :Scott)
# Plot
l = @layout [a{0.48w} _ b{0.48w}; c{0.48w} _ d{0.48w}]
plot(b_plt, c_plt, B_plt, C_plt, layout  = l, size = (1200,800), margin = 40Plots.mm)
savefig("figures/interpret/scott_example.pdf")

# Special case to show effect of C on H can be positive
param_sco_special = SilerParam(b = 5.0, B = 5.0, c = 0.01, C = 10.0, d = 0.)
Hgrad.([param_sco_special], 0:10, [:C], spec = :Scott)
sco_df_special = init_illus(param_sco_special, :Scott, ages = 0:40)
base_plt(sco_df_special, param_sco_special)
sco_df_special, C_plt = param_change_plt(sco_df_special, param_sco_special,
    :C; LE_change = 1, spec = :Scott)
C_plt
savefig("figures/interpret/scott_special_case.pdf")



"""
Bergeron specification example
"""
# Compute parameter changes
ber_df_1, b_plt = param_change_plt(ber_df_1[1:31,:], param_ber_1, :b; LE_change = 5, spec = :Bergeron)
ber_df_1, B_plt = param_change_plt(ber_df_1, param_ber_1, :B; LE_change = 5, spec = :Bergeron)
ber_df_2, c_plt = param_change_plt(ber_df_2, param_ber_2, :c; LE_change = 5, spec = :Bergeron)
ber_df_2, C_plt = param_change_plt(ber_df_2, param_ber_2, :C; LE_change = 5, spec = :Bergeron)
# Plot
l = @layout [a{0.48w} _ b{0.48w}; c{0.48w} _ d{0.48w}]
plot(b_plt, c_plt, B_plt, C_plt, layout  = l, size = (1200,800), margin = 40Plots.mm)
savefig("figures/interpret/bergeron_example.pdf")


# Special case to show effect of C on H can be positive
param_ber_special = SilerParam(b = 5.0, B = 0.00, c = 0.1, C = 75, d = 0.01)
plot(siler.([param_ber_special], 0:110, spec = :Bergeron), label = "Bergeron")
Hgrad.([param_ber_special], 0:10, [:C], spec = :Bergeron)
ber_df_special = init_illus(param_ber_special, :Bergeron)
base_plt(ber_df_special, param_ber_special)
ber_df_special, C_plt = param_change_plt(ber_df_special, param_ber_special,
    :C; LE_change = 5, spec = :Bergeron)
C_plt
savefig("figures/interpret/bergeron_special_case.pdf")
