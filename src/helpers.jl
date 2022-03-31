"""
Helper functions to analyse Turing results
"""


"""
Break matrix into chunks (array of arrays)
"""
chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]


"""
Get a positive definite covariance matrix
"""
quad_form_diag(M, v) = Symmetric((v .* v') .* (M .+ M') ./ 2)



"""
Function to extract summary statistics from time series variables
"""
function summarise_stats(df_sum, stats_post, years; parname = :unknown)

    rename!(df_sum, string.(years))
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
    # Add some convergence statistics
    df_out = hcat(df_ests, stats_post[:,[:mcse, :ess, :rhat, :ess_per_sec]])

    return(df_out)
end



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
Function that creates summary dataframes for each Siler parameter
    spec defines which Siler specification we want (Colchero, Scott, Bergeron or Standard)
    model_vers defines which dynamic model (indep, justrw, i2drift, firstdiff)
"""
function extract_variables(chain_in, years::Vector{Int64}; log_pars = false,
        spec = :Colchero, model_vers = :indep)

    # Siler parameter names
    if log_pars
        Bs = Symbol.("lB[".*string.(1:length(years)).*"]")
        bs = Symbol.("lb[".*string.(1:length(years)).*"]")
        Cs = Symbol.("lC[".*string.(1:length(years)).*"]")
        cs = Symbol.("lc[".*string.(1:length(years)).*"]")
        ds = Symbol.("ld[".*string.(1:length(years)).*"]")
        σs = Symbol.("lσ[".*string.(1:length(years)).*"]")
    else
        Bs = Symbol.("B[".*string.(1:length(years)).*"]")
        bs = Symbol.("b[".*string.(1:length(years)).*"]")
        Cs = Symbol.("C[".*string.(1:length(years)).*"]")
        cs = Symbol.("c[".*string.(1:length(years)).*"]")
        ds = Symbol.("d[".*string.(1:length(years)).*"]")
        σs = Symbol.("σ[".*string.(1:length(years)).*"]")
    end
    # Convert sampled chains to dataframe
    df_post = DataFrame(chain_in)
    # If parameters are logged, then convert now
    if log_pars
        df_post[:, Bs] = exp.(df_post[:, Bs])
        df_post[:, bs] = exp.(df_post[:, bs])
        df_post[:, Cs] = exp.(df_post[:, Cs])
        df_post[:, cs] = exp.(df_post[:, cs])
        df_post[:, ds] = exp.(df_post[:, ds])
        df_post[:, σs] = exp.(df_post[:, σs])
    end
    if spec == :Colchero
        df_post[:, Bs] = df_post[:, Bs]
        df_post[:, Cs] = df_post[:, Cs]
    elseif spec == :Scott
        df_post[:, Bs] = Matrix(df_post[:, Bs])./Matrix(df_post[:, bs])
        df_post[:, Cs] = Matrix(df_post[:, Cs])./Matrix(df_post[:, cs])
    elseif spec == :Bergeron
        df_post[:, Bs] = exp.(.- Matrix(df_post[:, Bs]))
        df_post[:, Cs] = (Matrix(df_post[:, Cs]) .+ log.(Matrix(df_post[:, cs])))./Matrix(df_post[:, cs])
    elseif spec == :Standard
        df_post[:, Bs] = exp.(.- Matrix(df_post[:, Bs]))
        df_post[:, Cs] = exp.(.- Matrix(df_post[:, Cs]))
    end


    stats = DataFrame(summarystats(chain_in))
    B_ests = summarise_stats(df_post[:,Bs], stats[in.(stats.parameters, [Bs]),:],
        years, parname = :B)
    b_ests = summarise_stats(df_post[:,bs], stats[in.(stats.parameters, [bs]),:],
        years, parname = :b)
    C_ests = summarise_stats(df_post[:,Cs], stats[in.(stats.parameters, [Cs]),:],
        years, parname = :C)
    c_ests = summarise_stats(df_post[:,cs], stats[in.(stats.parameters, [cs]),:],
        years, parname = :c)
    d_ests = summarise_stats(df_post[:,ds], stats[in.(stats.parameters, [ds]),:],
        years, parname = :d)
    σ_ests = summarise_stats(df_post[:,σs], stats[in.(stats.parameters, [σs]),:],
        years, parname = :σ)

    par_ests = vcat(B_ests, b_ests, C_ests, c_ests, d_ests, σ_ests)

    # If we have a dynamic model with various time series parameters, add these
    if model_vers == :justrw
        σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
        σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
            Int.(1:6), parname = :σ_par)
        σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
        σ_pars_ests.year .= 0
        # Add to dataframe
        par_ests = vcat(par_ests, σ_pars_ests)

    elseif model_vers == :firstdiff
        # Parameter time series variance
        σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
        σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
            Int.(1:6), parname = :σ_par)
        σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
        σ_pars_ests.year .= 0
        par_ests = vcat(par_ests, σ_pars_ests)
        # Parameter time series constant
        α_par_names = Symbol.("α_pars[".*string.(1:6).*"]")
        α_pars_ests = summarise_stats(df_post[:,α_par_names], stats[in.(stats.parameters, [α_par_names]),:],
            Int.(1:6), parname = :α_par)
        α_pars_ests.parameter = [:α_B, :α_b, :α_C, :α_c, :α_d, :α_σ]
        α_pars_ests.year .= 0
        par_ests = vcat(par_ests, α_pars_ests)
        # Parameter time series AR(1) coefficient
        β_par_names = Symbol.("β_pars[".*string.(1:6).*"]")
        β_pars_ests = summarise_stats(df_post[:,β_par_names], stats[in.(stats.parameters, [β_par_names]),:],
            Int.(1:6), parname = :β_par)
        β_pars_ests.parameter = [:β_B, :β_b, :β_C, :β_c, :β_d, :β_σ]
        β_pars_ests.year .= 0
        # Add to dataframe
        par_ests = vcat(par_ests, β_pars_ests)

    elseif model_vers == :i2drift
        σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
        σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
            Int.(1:6), parname = :σ_par)
        σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
        σ_pars_ests.year .= 0

        par_ests = vcat(par_ests, σ_pars_ests)

        σ_α_par_names = Symbol.("σ_α".*["B", "b", "C", "c", "d", "σ"])
        σ_α_pars_ests = summarise_stats(df_post[:,(σ_α_par_names)], stats[in.(stats.parameters, [σ_α_par_names]),:],
            Int.(1:6), parname = :σ_α_par)
        σ_α_pars_ests.parameter = [:σ_αB, :σ_αb, :σ_αC, :σ_αc, :σ_αd, :σ_ασ]
        σ_α_pars_ests.year .= 0

        par_ests = vcat(par_ests, σ_α_pars_ests)

        α_Bs = Symbol.("α_B[".*string.(1:length(years)-1).*"]")
        α_bs = Symbol.("α_b[".*string.(1:length(years)-1).*"]")
        α_Cs = Symbol.("α_C[".*string.(1:length(years)-1).*"]")
        α_cs = Symbol.("α_c[".*string.(1:length(years)-1).*"]")
        α_ds = Symbol.("α_d[".*string.(1:length(years)-1).*"]")
        α_σs = Symbol.("α_σ[".*string.(1:length(years)-1).*"]")

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

    end

    return par_ests

end





"""
Function to create a decomposed df from parameter estimates
"""
function create_decomp(parests_df; spec = :Colchero, eval_age = 0)

    # Convert from long to wide
    decomp_df = parests_df[.!occursin.("_",string.(parests_df.parameter)),
        ["year", "parameter", "median"]]
    decomp_df = unstack(decomp_df, :year, :parameter, :median)
    decomp_df = hcat(DataFrame(year = decomp_df[:,1]),Float64.(decomp_df[:,2:end]))
    # Make suer we're ordered chronologically
    decomp_df = sort!(decomp_df, :year)

    decomp_vars = ["LE_mod", "H_mod", "LE_b", "LE_B", "LE_c", "LE_C", "LE_d",
        "H_b", "H_B", "H_c", "H_C", "H_d", "Δb", "ΔB", "Δc", "ΔC", "Δd", "ΔLE_mod", "ΔH_mod"]
    decomp_df = hcat(decomp_df,DataFrame(NaN.*zeros(nrow(decomp_df), length(decomp_vars)), decomp_vars))

    # Get LE, H, changes and derivatives
    for ii in 1:nrow(decomp_df)

        # Identify parameters for that year
        row = NamedTuple(decomp_df[ii,:])
        params = SilerParam(b = row.b, B = row.B, c = row.c, C = row.C, d = row.d)
        # Compute model-impled life expectancy and inequality
        S_mod = siler_S.([params], [eval_age], eval_age:200, spec = spec)
        decomp_df.LE_mod[ii] = LE(params, eval_age, spec = spec)
        decomp_df.H_mod[ii] = H(params, eval_age, spec = spec)
        # Compute LE gradient for each parameter
        decomp_df.LE_B[ii] = LEgrad(params, eval_age, :B, spec = spec)
        decomp_df.LE_b[ii] = LEgrad(params, eval_age, :b, spec = spec)
        decomp_df.LE_C[ii] = LEgrad(params, eval_age, :C, spec = spec)
        decomp_df.LE_c[ii] = LEgrad(params, eval_age, :c, spec = spec)
        decomp_df.LE_d[ii] = LEgrad(params, eval_age, :d, spec = spec)
        # Compute H gradient for each parameter
        decomp_df.H_B[ii] = Hgrad(params, eval_age, :B, spec = spec)
        decomp_df.H_b[ii] = Hgrad(params, eval_age, :b, spec = spec)
        decomp_df.H_C[ii] = Hgrad(params, eval_age, :C, spec = spec)
        decomp_df.H_c[ii] = Hgrad(params, eval_age, :c, spec = spec)
        decomp_df.H_d[ii] = Hgrad(params, eval_age, :d, spec = spec)

        if ii > 1
            decomp_df.ΔB[ii] = decomp_df.B[ii] - decomp_df.B[ii-1]
            decomp_df.Δb[ii] = decomp_df.b[ii] - decomp_df.b[ii-1]
            decomp_df.ΔC[ii] = decomp_df.C[ii] - decomp_df.C[ii-1]
            decomp_df.Δc[ii] = decomp_df.c[ii] - decomp_df.c[ii-1]
            decomp_df.Δd[ii] = decomp_df.d[ii] - decomp_df.d[ii-1]
            decomp_df.ΔLE_mod[ii] = decomp_df.LE_mod[ii] - decomp_df.LE_mod[ii-1]
            decomp_df.ΔH_mod[ii] = decomp_df.H_mod[ii] - decomp_df.H_mod[ii-1]
        end

    end
    return decomp_df

end




"""
A function that extracts a vector of SilerParam structures from rows of a DataFrame
    of samples from a posterior distribution
"""
function particles2params(df_particle, Tptt; log_pars = true)
     B = Symbol("lB["*string(Tptt)*"]")
     b = Symbol("lb["*string(Tptt)*"]")
     C = Symbol("lC["*string(Tptt)*"]")
     c = Symbol("lc["*string(Tptt)*"]")
     d = Symbol("ld["*string(Tptt)*"]")

    params = repeat([SilerParam()], nrow(df_particle))
    if log_pars
        for ii in 1:nrow(df_particle)
            params[ii] = SilerParam(b = exp(df_particle[ii,b]), B = exp(df_particle[ii,B]),
                c = exp(df_particle[ii,c]), C = exp(df_particle[ii,C]),
                d = exp(df_particle[ii,d]))
        end
    else
        for ii in 1:nrow(df_particle)
            params[ii] = SilerParam(b = df_particle[ii,b], B = df_particle[ii,B],
                c = df_particle[ii,c], C = df_particle[ii,C], d = df_particle[ii,d])
        end
    end
    return params
end



## For each particle, simulate one step forward
"""
Function to compute predictive distribution nahead periods ahead with ndraws draws for
    the path of future shocks
"""
function compute_forecasts(df_post, nahead, ndraws, years; spec = :Colchero)

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

    # Extend dataframe to account for forward simulations
    df_pred = repeat(df_post, inner = ndraws)
    insertcols!(df_pred, 1, :particle => repeat([NaN], nrow(df_pred)) )
    df_pred[:,:particle] .= df_pred.iteration .+ df_pred.chain./10
    # Add extra columns for predicted values
    pred_vars = vcat(α_B_f, α_b_f, α_C_f, α_c_f, α_d_f, α_σ_f, B_f, b_f, C_f, c_f, d_f, σ_f)
    df_pred = hcat(df_pred,DataFrame(NaN.*zeros(nrow(df_pred), length(pred_vars)), pred_vars))
    # Split into groups for more efficient simulation
    df_grpd = groupby(df_pred, :particle)

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
            df_particle[:,σs[Tptt]] = df_particle[:,α_σs[Tptt-1]] + df_particle[:,σs[Tptt-1]] +
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




"""
Function that creates summary dataframes for current and forecast Siler parameters
    Need to specify which years are past and which are future

"""
function extract_forecast_variables(df_pred, past_years::Vector{Int64}, fut_years;
        log_pars = true, spec = :Colchero, model_vers = :i2drift)


    # Work on deep copy version as we might transform to account for logs
    df_in = deepcopy(df_pred)
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
        # If parameters are logged, then convert now
        df_in[:, Bs] = exp.(df_in[:, Bs])
        df_in[:, bs] = exp.(df_in[:, bs])
        df_in[:, Cs] = exp.(df_in[:, Cs])
        df_in[:, cs] = exp.(df_in[:, cs])
        df_in[:, ds] = exp.(df_in[:, ds])
        df_in[:, σs] = exp.(df_in[:, σs])
    else
        Bs = Symbol.("B[".*string.(1:length(all_years)).*"]")
        bs = Symbol.("b[".*string.(1:length(all_years)).*"]")
        Cs = Symbol.("C[".*string.(1:length(all_years)).*"]")
        cs = Symbol.("c[".*string.(1:length(all_years)).*"]")
        ds = Symbol.("d[".*string.(1:length(all_years)).*"]")
        σs = Symbol.("σ[".*string.(1:length(all_years)).*"]")
    end

    if spec == :Colchero
        df_in[:, Bs] = df_in[:, Bs]
        df_in[:, Cs] = df_in[:, Cs]
    elseif spec == :Scott
        df_in[:, Bs] = Matrix(df_in[:, Bs])./Matrix(df_in[:, bs])
        df_in[:, Cs] = Matrix(df_in[:, Cs])./Matrix(df_in[:, cs])
    elseif spec == :Bergeron
        df_in[:, Bs] = exp.(.- Matrix(df_in[:, Bs]))
        df_in[:, Cs] = (Matrix(df_in[:, Cs]) .+ log.(Matrix(df_in[:, cs])))./Matrix(df_in[:, cs])
    elseif spec == :Standard
        df_in[:, Bs] = exp.(.- Matrix(df_in[:, Bs]))
        df_in[:, Cs] = exp.(.- Matrix(df_in[:, Cs]))
    end


    B_ests = summarise_forecasts(df_in[:,Bs], all_years, parname = :B)
    b_ests = summarise_forecasts(df_in[:,bs], all_years, parname = :b)
    C_ests = summarise_forecasts(df_in[:,Cs], all_years, parname = :C)
    c_ests = summarise_forecasts(df_in[:,cs], all_years, parname = :c)
    d_ests = summarise_forecasts(df_in[:,ds], all_years, parname = :d)
    σ_ests = summarise_forecasts(df_in[:,σs], all_years, parname = :σ)

    par_ests = vcat(B_ests, b_ests, C_ests, c_ests, d_ests, σ_ests)

    # If we have a dynamic model with various time series parameters, add these
    if model_vers == :justrw

    elseif model_vers == :firstdiff

    elseif model_vers == :i2drift

        σ_par_names = Symbol.("σ_pars[".*string.(1:6).*"]")
        σ_pars_ests = summarise_forecasts(df_in[:,σ_par_names],Int.(1:6), parname = :σ_par)
        σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
        σ_pars_ests.year .= 0
        par_ests = vcat(par_ests, σ_pars_ests)

        σ_α_par_names = Symbol.("σ_α".*["B", "b", "C", "c", "d", "σ"])
        σ_α_pars_ests = summarise_forecasts(df_in[:,(σ_α_par_names)], Int.(1:6), parname = :σ_α_par)
        σ_α_pars_ests.parameter = [:σ_αB, :σ_αb, :σ_αC, :σ_αc, :σ_αd, :σ_ασ]
        σ_α_pars_ests.year .= 0
        par_ests = vcat(par_ests, σ_α_pars_ests)

        α_Bs = Symbol.("α_B[".*string.(1:length(all_years)-1).*"]")
        α_bs = Symbol.("α_b[".*string.(1:length(all_years)-1).*"]")
        α_Cs = Symbol.("α_C[".*string.(1:length(all_years)-1).*"]")
        α_cs = Symbol.("α_c[".*string.(1:length(all_years)-1).*"]")
        α_ds = Symbol.("α_d[".*string.(1:length(all_years)-1).*"]")
        α_σs = Symbol.("α_σ[".*string.(1:length(all_years)-1).*"]")

        α_B_ests = summarise_forecasts(df_in[:,α_Bs], all_years[2:end], parname = :α_B)
        α_b_ests = summarise_forecasts(df_in[:,α_bs], all_years[2:end], parname = :α_b)
        α_C_ests = summarise_forecasts(df_in[:,α_Cs], all_years[2:end], parname = :α_C)
        α_c_ests = summarise_forecasts(df_in[:,α_cs], all_years[2:end], parname = :α_c)
        α_d_ests = summarise_forecasts(df_in[:,α_ds], all_years[2:end], parname = :α_d)
        α_σ_ests = summarise_forecasts(df_in[:,α_σs], all_years[2:end], parname = :α_σ)

        par_ests = vcat(par_ests, α_B_ests, α_b_ests, α_C_ests, α_c_ests, α_d_ests, α_σ_ests)

    end

    # Extrac the model implied forecasts of LE and H
    LEs = Symbol.("LE[".*string.(1:length(all_years)).*"]")
    Hs = Symbol.("H[".*string.(1:length(all_years)).*"]")

    LE_ests = summarise_forecasts(df_in[:,LEs], all_years, parname = :LE)
    H_ests = summarise_forecasts(df_in[:,Hs], all_years, parname = :H)

    par_ests = vcat(par_ests, LE_ests, H_ests)


    insertcols!(par_ests, 2, :forecast => repeat([0], nrow(par_ests)) )
    par_ests.forecast[in.(par_ests.year, [vcat(fut_years, past_years[end])])] .= 1

    return par_ests

end
