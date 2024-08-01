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
Import some data
"""
all_df = CSV.read("data/clean/all_lifetab_5y.csv", DataFrame, ntasks = 1);
bp_df = all_df[(all_df.best_practice .== 1) .& (all_df.year .>= 1900), :];
sort!(bp_df, [:year, :age]);


"""
Select parameters
"""
## Parameters to look at infant mortality
param_col_1 = SilerParam(b = 1.0, B = 1.25, c = 0.15, C = 10.0, d = 0.001)
param_sco_1 = SilerParam(b = param_col_1.b, B = param_col_1.B/param_col_1.b,
    c = param_col_1.c, C = param_col_1.C/param_col_1.c, d = param_col_1.d)
param_ber_1 = SilerParam(b = param_col_1.b, B= (param_col_1.B + log(param_col_1.b))/param_col_1.b,
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
param_ber_2 = SilerParam(b = param_col_2.b, B= (param_col_2.B + log(param_col_2.b))/param_col_2.b,
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

# Figure for the paper
ber_df_c, c_plt = param_change_plt(ber_df_2, param_ber_2, :c; LE_change = 5, spec = :Bergeron)
ber_df_C, C_plt = param_change_plt(ber_df_2, param_ber_2, :C; LE_change = 5, spec = :Bergeron)
plot(ber_df_c.age, ber_df_c.S_base, label = "Baseline", ylabel = "Survival Rate", xlabel = "Age")
plot!(ber_df_c.age, ber_df_c.S_cchange, linestyle = :dash, label = "Change to c")
plot!(ber_df_C.age, ber_df_C.S_Cchange, linestyle = :dash, label = "Change to C")
plot!(legend = :bottomleft, size = (400,300))
savefig("figures/interpret/survival_c_example.pdf")

ber_df_b, b_plt = param_change_plt(ber_df_1, param_ber_1, :b; LE_change = 5, spec = :Bergeron)
ber_df_B, B_plt = param_change_plt(ber_df_1, param_ber_1, :B; LE_change = 5, spec = :Bergeron)
plot(ber_df_b.age, ber_df_b.S_base, label = "Baseline", ylabel = "Survival Rate", xlabel = "Age")
plot!(ber_df_b.age, ber_df_b.S_bchange, linestyle = :dash, label = "Change to b")
plot!(ber_df_B.age, ber_df_B.S_Bchange, linestyle = :dash, label = "Change to B")
plot!(legend = :bottomleft, size = (400,300))
savefig("figures/interpret/survival_b_example.pdf")

plot(ber_df_b.age, ber_df_b.μ_base, label = "Baseline", ylabel = "Survival Rate", xlabel = "Age")
plot!(ber_df_b.age, ber_df_b.μ_bchange, linestyle = :dash, label = "Change to b")
plot!(ber_df_B.age, ber_df_B.μ_Bchange, linestyle = :dash, label = "Change to B")
plot!(legend = :topright, size = (400,300))
savefig("figures/interpret/mortality_b_example.pdf")



# Special case to show effect of C on H can be positive
param_ber_special = SilerParam(b = 1.0, B = 5.0, c = 0.1, C = 75, d = 0.01)
plot(siler.([param_ber_special], 0:110, spec = :Bergeron), label = "Bergeron")
Hgrad.([param_ber_special], 0:10, [:C], spec = :Bergeron)
ber_df_special = init_illus(param_ber_special, :Bergeron)
base_plt(ber_df_special, param_ber_special)
ber_df_special, C_plt = param_change_plt(ber_df_special, param_ber_special,
    :C; LE_change = 5, spec = :Bergeron)
C_plt
savefig("figures/interpret/bergeron_special_case.pdf")



"""
Heligman Pollard example
"""
param_hp = HPParam(A = 0.0005, B = 0.04, C = 0.12, D = 0.05, E = 6.0,
    F = 25.0, G = 0.00009, H = 1.08, σ = 0.0)

mort_rates = HeligmanPollard.([param_hp], 0:100)
plot(log.(mort_rates))

sum(cumprod(1 .- mort_rates))

ages = 0:109
first_part = param_hp.A.^((ages .+ param_hp.B).^param_hp.C) 
second_part = param_hp.D.*exp.(.-param_hp.E.*(log.(ages) .- log.(param_hp.F)).^2) 
third_part = param_hp.G.*param_hp.H.^(ages) 

plot(log.(first_part), label = "Infant")
plot!(log.(second_part), label = "Accident")
plot!(log.(third_part), label = "Later-life")
plot!(log.(first_part + second_part + third_part), label = "Total")


# Extract BP mortality data
country_df = bp_df;
country_df[country_df[:,:mx_f] .== 0.0,:mx_f] .=  minimum(country_df[country_df[:,:mx_f] .> 0.0,:mx_f]);
country_m_data = chunk(country_df.mx_f, maximum(country_df.age)+1);
country_lm_data = [log.(m_dist) for m_dist in country_m_data];
plot(country_lm_data, legend = false)
country_ages = Int64.(0:maximum(country_df.age))
country_years = unique(country_df.year)
T = length(country_lm_data)
@assert length(country_m_data)==length(country_years) "number of years doesn't match length of m_data"
@assert length(country_m_data[1])==length(country_ages) "number of ages doesn't match length of m_data[1]"

# Standard static Siler model
niters = 100
nchains = 1
tt = 24
prior_hp = sample(log_hp_static(country_lm_data[tt], country_ages), Prior(), 50)


param_hp = HPParam(A = 0.0005, B = 0.04, C = 0.12, D = 0.05, E = 6.0,
    F = 25.0, G = 0.00009, H = 1.08, σ = 0.0)

mort_rates = HeligmanPollard.([param_hp], 0:100)
plot(log.(mort_rates))


chain_static = sample(log_hp_static(country_lm_data[tt], country_ages), NUTS(0.65), MCMCThreads(),
    niters, nchains)
display(chain_static)





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


# Figure for the paper
col_df_c, c_plt = param_change_plt(col_df_2, param_col_2, :c; LE_change = 5, spec = :Colchero)
col_df_C, C_plt = param_change_plt(col_df_2, param_col_2, :C; LE_change = 5, spec = :Colchero)
plot(col_df_c.age, col_df_c.S_base, label = "Baseline", ylabel = "Survival Rate", xlabel = "Age")
plot!(col_df_c.age, col_df_c.S_cchange, linestyle = :dash, label = "Change to c")
plot!(col_df_C.age, col_df_C.S_Cchange, linestyle = :dash, label = "Change to C")
plot!(legend = :bottomleft, size = (400,300))
savefig("figures/interpret/survival_c_example_col.pdf")


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
param_sco_special = SilerParam(b = 5.0, B = 5.0, c = 0.1, C = 100.0, d = 0.01)
Hgrad.([param_sco_special], 0:10, [:C], spec = :Scott)
sco_df_special = init_illus(param_sco_special, :Scott)
base_plt(sco_df_special, param_sco_special)
sco_df_special, C_plt = param_change_plt(sco_df_special, param_sco_special,
    :C; LE_change = 5, spec = :Scott)
C_plt
savefig("figures/interpret/scott_special_case.pdf")
