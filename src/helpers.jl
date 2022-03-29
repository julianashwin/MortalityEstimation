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
