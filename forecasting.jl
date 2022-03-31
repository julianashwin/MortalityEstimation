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
folder = "figures/best_practice/"
model = "i2"
nahead = 1

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



"""
Identify the variable names
"""
## Siler parameters
Bs = Symbol.("lB[".*string.(1:length(years)+nahead).*"]")
bs = Symbol.("lb[".*string.(1:length(years)+nahead).*"]")
Cs = Symbol.("lC[".*string.(1:length(years)+nahead).*"]")
cs = Symbol.("lc[".*string.(1:length(years)+nahead).*"]")
ds = Symbol.("ld[".*string.(1:length(years)+nahead).*"]")
σs = Symbol.("lσ[".*string.(1:length(years)+nahead).*"]")
## Drift terms (remeber lB[tt] = α_B[tt-1] + lB[tt-1] = shock)
α_Bs = Symbol.("α_B[".*string.(1:length(years)-1+nahead).*"]")
α_bs = Symbol.("α_b[".*string.(1:length(years)-1+nahead).*"]")
α_Cs = Symbol.("α_C[".*string.(1:length(years)-1+nahead).*"]")
α_cs = Symbol.("α_c[".*string.(1:length(years)-1+nahead).*"]")
α_ds = Symbol.("α_d[".*string.(1:length(years)-1+nahead).*"]")
α_σs = Symbol.("α_σ[".*string.(1:length(years)-1+nahead).*"]")
## Variance terms
σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
plot(chain_post[σ_par_names])
σ_α_par_names = Symbol.("σ_α".*["B", "b", "C", "c", "d", "σ"])
plot(chain_post[σ_α_par_names])

## The variables to be forecast
# Future parameters
B_f = Symbol.("lB[".*string.(length(years)+1:length(years)+nahead).*"]")
b_f = Symbol.("lb[".*string.(length(years)+1:length(years)+nahead).*"]")
C_f = Symbol.("lC[".*string.(length(years)+1:length(years)+nahead).*"]")
c_f = Symbol.("lc[".*string.(length(years)+1:length(years)+nahead).*"]")
d_f = Symbol.("ld[".*string.(length(years)+1:length(years)+nahead).*"]")
σ_f = Symbol.("lσ[".*string.(length(years)+1:length(years)+nahead).*"]")
# Future drift terms
α_B_f = Symbol.("α_B[".*string.(length(years):length(years)-1+nahead).*"]")
α_b_f = Symbol.("α_b[".*string.(length(years):length(years)-1+nahead).*"]")
α_C_f = Symbol.("α_C[".*string.(length(years):length(years)-1+nahead).*"]")
α_c_f = Symbol.("α_c[".*string.(length(years):length(years)-1+nahead).*"]")
α_d_f = Symbol.("α_d[".*string.(length(years):length(years)-1+nahead).*"]")
α_σ_f = Symbol.("α_σ[".*string.(length(years):length(years)-1+nahead).*"]")

## The model-implied variables
LEs = Symbol.("LE[".*string.(1:length(years)+nahead).*"]")
Hs = Symbol.("H[".*string.(1:length(years)+nahead).*"]")


## Compute the model implied LE and H for the in-sample periods
# Add extra columns for model impled LE and H
impl_vars = vcat(LEs, Hs)
df_post = hcat(df_post,DataFrame(NaN.*zeros(nrow(df_post), length(impl_vars)), impl_vars))
# Loop through each period
prog = Progress(length(years), desc = "Calculating model implied LE and H: ")
for ii in 1:length(years)
    params = particles2params(df_post, ii, Bs, bs, Cs, cs, ds; log_pars = true)
    df_post[:, LEs[ii]] = LE.(params, [0.0], spec = :Colchero)
    df_post[:, Hs[ii]] = H.(params, [0.0], spec = :Colchero)
    next!(prog)
end



"""
Create larger DataFrame for predictive distribution
"""
## Set up some infrastructure
# Define number of draws to approximate each future shock
ndraws = 10
# 
# Split into groups for more efficient simulation
df_grpd = groupby(df_pred, :particle)



## For each particle, simulate one step forward
function compute_forecasts(df_post, nahead, ndraws, years; spec = :Colchero)

    # Extend dataframe to account for forward simulations
    df_pred = repeat(df_post, inner = ndraws*nahead)
    insertcols!(df_pred, 1, :particle => repeat([NaN], nrow(df_pred)) )
    df_pred[:,:particle] .= df_pred.iteration .+ df_pred.chain./10
    # Add extra columns for predicted values
    pred_vars = vcat(α_B_f, α_b_f, α_C_f, α_c_f, α_d_f, α_σ_f, B_f, b_f, C_f, c_f, d_f, σ_f)
    df_pred = hcat(df_pred,DataFrame(NaN.*zeros(nrow(df_pred), length(pred_vars)), pred_vars))

    # Siler parameters
    Bs = Symbol.("lB[".*string.(1:length(years)+nahead).*"]")
    bs = Symbol.("lb[".*string.(1:length(years)+nahead).*"]")
    Cs = Symbol.("lC[".*string.(1:length(years)+nahead).*"]")
    cs = Symbol.("lc[".*string.(1:length(years)+nahead).*"]")
    ds = Symbol.("ld[".*string.(1:length(years)+nahead).*"]")
    σs = Symbol.("lσ[".*string.(1:length(years)+nahead).*"]")
    # Drift terms (remeber lB[tt] = α_B[tt-1] + lB[tt-1] = shock)
    α_Bs = Symbol.("α_B[".*string.(1:length(years)-1+nahead).*"]")
    α_bs = Symbol.("α_b[".*string.(1:length(years)-1+nahead).*"]")
    α_Cs = Symbol.("α_C[".*string.(1:length(years)-1+nahead).*"]")
    α_cs = Symbol.("α_c[".*string.(1:length(years)-1+nahead).*"]")
    α_ds = Symbol.("α_d[".*string.(1:length(years)-1+nahead).*"]")
    α_σs = Symbol.("α_σ[".*string.(1:length(years)-1+nahead).*"]")
    # Future parameters
    B_f = Symbol.("lB[".*string.(length(years)+1:length(years)+nahead).*"]")
    b_f = Symbol.("lb[".*string.(length(years)+1:length(years)+nahead).*"]")
    C_f = Symbol.("lC[".*string.(length(years)+1:length(years)+nahead).*"]")
    c_f = Symbol.("lc[".*string.(length(years)+1:length(years)+nahead).*"]")
    d_f = Symbol.("ld[".*string.(length(years)+1:length(years)+nahead).*"]")
    σ_f = Symbol.("lσ[".*string.(length(years)+1:length(years)+nahead).*"]")
    # Future drift terms
    α_B_f = Symbol.("α_B[".*string.(length(years):length(years)-1+nahead).*"]")
    α_b_f = Symbol.("α_b[".*string.(length(years):length(years)-1+nahead).*"]")
    α_C_f = Symbol.("α_C[".*string.(length(years):length(years)-1+nahead).*"]")
    α_c_f = Symbol.("α_c[".*string.(length(years):length(years)-1+nahead).*"]")
    α_d_f = Symbol.("α_d[".*string.(length(years):length(years)-1+nahead).*"]")
    α_σ_f = Symbol.("α_σ[".*string.(length(years):length(years)-1+nahead).*"]")
    # The model-implied variables
    LEs = Symbol.("LE[".*string.(1:length(years)+nahead).*"]")
    Hs = Symbol.("H[".*string.(1:length(years)+nahead).*"]")

    # Loop through years ahead for forecasts
    T = length(years)
    for tt in 1:nahead
        Tptt = T + tt

        # Loop through each particle
        prog = Progress(length(df_grpd), desc = "Simulating "*string(tt)*" periods ahead: ")
        for gg in 1:length(df_grpd)
            df_particle = df_grpd[gg]
            # Parameter variance terms for this particle
            σ_ϵB = df_particle[1,Symbol("σ_pars[1]")]
            σ_ϵb = df_particle[1,Symbol("σ_pars[2]")]
            σ_ϵC = df_particle[1,Symbol("σ_pars[3]")]
            σ_ϵc = df_particle[1,Symbol("σ_pars[4]")]
            σ_ϵd = df_particle[1,Symbol("σ_pars[5]")]
            σ_ϵσ = df_particle[1,Symbol("σ_pars[6]")]
            # Drift variance terms for this particle
            σ_ξB = df_particle[1,Symbol("σ_αB")]
            σ_ξb = df_particle[1,Symbol("σ_αb")]
            σ_ξC = df_particle[1,Symbol("σ_αC")]
            σ_ξc = df_particle[1,Symbol("σ_αc")]
            σ_ξd = df_particle[1,Symbol("σ_αd")]
            σ_ξσ = df_particle[1,Symbol("σ_ασ")]

            # Forecast drift terms (remeber that we need to go one further back for these)
            df_particle[:,α_Bs[Tptt-1]] = df_particle[:,α_Bs[Tptt-2]] + rand(Normal(0.0, σ_ξB),ndraws)
            df_particle[:,α_bs[Tptt-1]] = df_particle[:,α_bs[Tptt-2]] + rand(Normal(0.0, σ_ξb),ndraws)
            df_particle[:,α_Cs[Tptt-1]] = df_particle[:,α_Cs[Tptt-2]] + rand(Normal(0.0, σ_ξC),ndraws)
            df_particle[:,α_cs[Tptt-1]] = df_particle[:,α_cs[Tptt-2]] + rand(Normal(0.0, σ_ξc),ndraws)
            df_particle[:,α_ds[Tptt-1]] = df_particle[:,α_ds[Tptt-2]] + rand(Normal(0.0, σ_ξd),ndraws)
            df_particle[:,α_σs[Tptt-1]] = df_particle[:,α_σs[Tptt-2]] + rand(Normal(0.0, σ_ξσ),ndraws)

            # Forecast the parameters
            df_particle[:,Bs[Tptt]] = df_particle[:,α_Bs[Tptt-1]] + df_particle[:,Bs[Tptt-1]] +
                rand(Normal(0.0, σ_ϵB),ndraws)
            df_particle[:,bs[Tptt]] = df_particle[:,α_bs[Tptt-1]] + df_particle[:,bs[Tptt-1]] +
                rand(Normal(0.0, σ_ϵb),ndraws)
            df_particle[:,Cs[Tptt]] = df_particle[:,α_Cs[Tptt-1]] + df_particle[:,Cs[Tptt-1]] +
                rand(Normal(0.0, σ_ϵC),ndraws)
            df_particle[:,cs[Tptt]] = df_particle[:,α_cs[Tptt-1]] + df_particle[:,cs[Tptt-1]] +
                rand(Normal(0.0, σ_ϵc),ndraws)
            df_particle[:,ds[Tptt]] = df_particle[:,α_ds[Tptt-1]] + df_particle[:,ds[Tptt-1]] +
                rand(Normal(0.0, σ_ϵd),ndraws)
            df_particle[:,σs[Tptt]] = df_particle[:,α_ds[Tptt-1]] + df_particle[:,ds[Tptt-1]] +
                rand(Normal(0.0, σ_ϵd),ndraws)

            # Compute the model implied LE and H overtime
            params = particles2params(df_particle, Tptt, Bs, bs, Cs, cs, ds; log_pars = true)
            df_particle[:,LEs[Tptt]] = LE.(params, [0.0], spec = spec)
            df_particle[:,Hs[Tptt]] = H.(params, [0.0], spec = spec)

            next!(prog)
        end
    end

    return df_pred
end


past_years = years
fut_years = Int.(maximum(years) .+ 5.0.*ones(nahead))
df_pred = extract_forecast_variables(df_pred, past_years, fut_years,
        log_pars = true, spec = :Colchero, model_vers = :i2drift)
