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


"""
Create larger DataFrame for predictive distribution
"""
## Set up some infrastructure
# Define number of draws to approximate each future shock
ndraws = 10
# Extend dataframe to account for forward simulations
df_pred = repeat(df_post, inner = ndraws*nahead)
insertcols!(df_pred, 1, :particle => repeat([NaN], nrow(df_pred)) )
df_pred[:,:particle] .= df_pred.iteration .+ df_pred.chain./10
# Add extra columns for predicted values
pred_vars = vcat(α_B_f, α_b_f, α_C_f, α_c_f, α_d_f, α_σ_f, B_f, b_f, C_f, c_f, d_f, σ_f)
df_pred = hcat(df_pred,DataFrame(NaN.*zeros(nrow(df_pred), length(pred_vars)), pred_vars))
# Add extra columns for model impled LE and H
impl_vars = vcat(LEs, Hs)
df_pred = hcat(df_pred,DataFrame(NaN.*zeros(nrow(df_pred), length(impl_vars)), impl_vars))
# Split into groups for more efficient simulation
df_grpd = groupby(df_pred, :particle)


## For each particle, simulate one step forward
T = length(years)
for tt in 1:nahead
    Tptt = T + tt

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
        params = SilerParam.(b = df_particle[:,bs[Tptt]], B = df_particle[:,Bs[Tptt]],
            c = df_particle[:,cs[Tptt]], C = df_particle[:,Cs[Tptt]], d = df_particle[:,ds[Tptt]])
        df_particle[:,LEs[Tptt]]

        next!(prog)
    end
end


function particles2params(df_particle, Tptt, spec = :Colchero)
    SilerParam(b = df_particle[:,Bs[Tptt]])

past_years = years
fut_years = Int.(maximum(years) .+ 5.0.*ones(nahead))


"""
Function to extract summary statistics from time series forecasts
"""
function summarise_forecasts(df_sum, all_years; parname = :unknown)

    rename!(df_sum, string.(all_years))
    df_ests = describe(df_sum)
    rename!(df_ests, Dict(:variable => :year))
    df_ests.year = parse.(Int64,string.(df_ests.year))
    insertcols!(df_ests, 1, :parameter => repeat([parname], nrow(df_ests)) )

    # Add some percentiles
    df_ests[:,:std] .= 0.0
    df_ests[:,:pc975] .= 0.0
    df_ests[:,:pc025] .= 0.0
    df_ests[:,:pc85] .= 0.0
    df_ests[:,:pc15] .= 0.0
    df_ests[:,:pc75] .= 0.0
    df_ests[:,:pc25] .= 0.0
    for yy in 1:ncol(df_sum)
        df_ests.std[(yy)] = std(df_sum[:,yy])
        df_ests.pc975[(yy)] = percentile(df_sum[:,yy], 97.5)
        df_ests.pc025[(yy)] = percentile(df_sum[:,yy], 2.5)
        df_ests.pc85[(yy)] = percentile(df_sum[:,yy], 85)
        df_ests.pc15[(yy)] = percentile(df_sum[:,yy], 15)
        df_ests.pc75[(yy)] = percentile(df_sum[:,yy], 75)
        df_ests.pc25[(yy)] = percentile(df_sum[:,yy], 25)
    end

    return(df_ests)
end



"""
Function that creates summary dataframes for current and forecast Siler parameters
    Need to specify which years are past and which are future

"""
function extract_forecast_variables(df_pred, past_years::Vector{Int64}, fut_years;
        log_pars = false, spec = :Colchero, model_vers = :i2drift)

    if model_vers != :i2drift
        throw("Only supported for i2drift")
    end
    all_years = vcat(past_years, fut_years)
    # Siler parameter names
    if log_pars
        Bs = Symbol.("lB[".*string.(1:length(all_years)).*"]")
        bs = Symbol.("lb[".*string.(1:length(all_years)).*"]")
        Cs = Symbol.("lC[".*string.(1:length(all_years)).*"]")
        cs = Symbol.("lc[".*string.(1:length(all_years)).*"]")
        ds = Symbol.("ld[".*string.(1:length(all_years)).*"]")
        σs = Symbol.("lσ[".*string.(1:length(all_years)).*"]")
    else
        Bs = Symbol.("B[".*string.(1:length(all_years)).*"]")
        bs = Symbol.("b[".*string.(1:length(all_years)).*"]")
        Cs = Symbol.("C[".*string.(1:length(all_years)).*"]")
        cs = Symbol.("c[".*string.(1:length(all_years)).*"]")
        ds = Symbol.("d[".*string.(1:length(all_years)).*"]")
        σs = Symbol.("σ[".*string.(1:length(all_years)).*"]")
    end
    # If parameters are logged, then convert now
    if log_pars
        df_pred[:, Bs] = exp.(df_pred[:, Bs])
        df_pred[:, bs] = exp.(df_pred[:, bs])
        df_pred[:, Cs] = exp.(df_pred[:, Cs])
        df_pred[:, cs] = exp.(df_pred[:, cs])
        df_pred[:, ds] = exp.(df_pred[:, ds])
        df_pred[:, σs] = exp.(df_pred[:, σs])
    end
    if spec == :Colchero
        df_pred[:, Bs] = df_pred[:, Bs]
        df_pred[:, Cs] = df_pred[:, Cs]
    elseif spec == :Scott
        df_pred[:, Bs] = Matrix(df_pred[:, Bs])./Matrix(df_pred[:, bs])
        df_pred[:, Cs] = Matrix(df_pred[:, Cs])./Matrix(df_pred[:, cs])
    elseif spec == :Bergeron
        df_pred[:, Bs] = exp.(.- Matrix(df_pred[:, Bs]))
        df_pred[:, Cs] = (Matrix(df_pred[:, Cs]) .+ log.(Matrix(df_pred[:, cs])))./Matrix(df_pred[:, cs])
    elseif spec == :Standard
        df_pred[:, Bs] = exp.(.- Matrix(df_pred[:, Bs]))
        df_pred[:, Cs] = exp.(.- Matrix(df_pred[:, Cs]))
    end


    B_ests = summarise_forecasts(df_pred[:,Bs], all_years, parname = :B)
    b_ests = summarise_forecasts(df_pred[:,bs], all_years, parname = :b)
    C_ests = summarise_forecasts(df_pred[:,Cs], all_years, parname = :C)
    c_ests = summarise_forecasts(df_pred[:,cs], all_years, parname = :c)
    d_ests = summarise_forecasts(df_pred[:,ds], all_years, parname = :d)
    σ_ests = summarise_forecasts(df_pred[:,σs], all_years, parname = :σ)

    par_ests = vcat(B_ests, b_ests, C_ests, c_ests, d_ests, σ_ests)

    # If we have a dynamic model with various time series parameters, add these
    if model_vers == :justrw

    elseif model_vers == :firstdiff

    elseif model_vers == :i2drift
        α_Bs = Symbol.("α_B[".*string.(1:length(all_years)-1).*"]")
        α_bs = Symbol.("α_b[".*string.(1:length(all_years)-1).*"]")
        α_Cs = Symbol.("α_C[".*string.(1:length(all_years)-1).*"]")
        α_cs = Symbol.("α_c[".*string.(1:length(all_years)-1).*"]")
        α_ds = Symbol.("α_d[".*string.(1:length(all_years)-1).*"]")
        α_σs = Symbol.("α_σ[".*string.(1:length(all_years)-1).*"]")

        α_B_ests = summarise_forecasts(df_pred[:,α_Bs], all_years[2:end], parname = :α_B)
        α_b_ests = summarise_forecasts(df_pred[:,α_bs], all_years[2:end], parname = :α_b)
        α_C_ests = summarise_forecasts(df_pred[:,α_Cs], all_years[2:end], parname = :α_C)
        α_c_ests = summarise_forecasts(df_pred[:,α_cs], all_years[2:end], parname = :α_c)
        α_d_ests = summarise_forecasts(df_pred[:,α_ds], all_years[2:end], parname = :α_d)
        α_σ_ests = summarise_forecasts(df_pred[:,α_σs], all_years[2:end], parname = :α_σ)

        par_ests = vcat(par_ests, α_B_ests, α_b_ests, α_C_ests, α_c_ests, α_d_ests, α_σ_ests)

    end

    insertcols!(par_ests, 2, :forecast => repeat([0], nrow(par_ests)) )
    par_ests.forecast[in.(par_ests.year, fut_years)] .= 1

    return par_ests

end

df_pred
