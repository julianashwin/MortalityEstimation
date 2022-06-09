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
        df_ests.std[(yy)] = std(df_sum[.!isnan.(df_sum[:,yy]),yy])
        df_ests.pc975[(yy)] = percentile(skipmissing(df_sum[.!isnan.(df_sum[:,yy]),yy]), 97.5)
        df_ests.pc025[(yy)] = percentile(df_sum[.!isnan.(df_sum[:,yy]),yy], 2.5)
        df_ests.pc85[(yy)] = percentile(df_sum[.!isnan.(df_sum[:,yy]),yy], 85)
        df_ests.pc15[(yy)] = percentile(df_sum[.!isnan.(df_sum[:,yy]),yy], 15)
        df_ests.pc75[(yy)] = percentile(df_sum[.!isnan.(df_sum[:,yy]),yy], 75)
        df_ests.pc25[(yy)] = percentile(df_sum[.!isnan.(df_sum[:,yy]),yy], 25)
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
        if model_vers == :cov
            Bs = Symbol.("lpars[".*string.(1:length(years)).*"][1]")
            bs = Symbol.("lpars[".*string.(1:length(years)).*"][2]")
            Cs = Symbol.("lpars[".*string.(1:length(years)).*"][3]")
            cs = Symbol.("lpars[".*string.(1:length(years)).*"][4]")
            ds = Symbol.("lpars[".*string.(1:length(years)).*"][5]")
            σs = Symbol.("lpars[".*string.(1:length(years)).*"][6]")
        else
            Bs = Symbol.("lB[".*string.(1:length(years)).*"]")
            bs = Symbol.("lb[".*string.(1:length(years)).*"]")
            Cs = Symbol.("lC[".*string.(1:length(years)).*"]")
            cs = Symbol.("lc[".*string.(1:length(years)).*"]")
            ds = Symbol.("ld[".*string.(1:length(years)).*"]")
            σs = Symbol.("lσ[".*string.(1:length(years)).*"]")
        end
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
        df_post[:, Bs] = Matrix(df_post[:, bs]).*Matrix(df_post[:, Bs]) .- log.(Matrix(df_post[:, bs]))
        df_post[:, Cs] = Matrix(df_post[:, cs]).*Matrix(df_post[:, Cs]) .- log.(Matrix(df_post[:, cs]))
    elseif spec == :Scott
        df_post[:, Bs] = Matrix(df_post[:, Bs]) .- log.(Matrix(df_post[:, bs]))./Matrix(df_post[:, bs])
        df_post[:, Cs] = Matrix(df_post[:, Cs]) .- log.(Matrix(df_post[:, cs]))./Matrix(df_post[:, cs])
    elseif spec == :Bergeron
        df_post[:, Bs] = df_post[:, Bs]
        df_post[:, Cs] = df_post[:, Cs]
    elseif spec == :Standard
        throw("Haven't worked out the new conversion for Standard yet")
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

    elseif (model_vers == :i2drift) | (model_vers == :cov)

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

        if model_vers == :cov
            ρ_s = [:ρ_bB, :ρ_cC, :ρ_cd, :ρ_Cd]
            ρ_ests = summarise_stats(df_post[:,ρ_s], stats[in.(stats.parameters, [ρ_s]),:],
                Int.(1:length(ρ_s)), parname = :ρ)
            ρ_ests.parameter = [:ρ_bB, :ρ_cC, :ρ_cd, :ρ_Cd]
            ρ_ests.year .= 0

            par_ests = vcat(par_ests, ρ_ests)
        end

    end

    return par_ests

end





"""
Function to create a decomposed df from parameter estimates
"""
function create_decomp(parests_df; spec = :Bergeron, eval_age = 0, forecasts = false)

    # Convert from long to wide
    decomp_df = parests_df[.!occursin.("_",string.(parests_df.parameter)),
        ["year", "parameter", "median"]]
    decomp_df = unstack(decomp_df, :year, :parameter, :median)
    decomp_df = hcat(DataFrame(year = decomp_df[:,1]),Float64.(decomp_df[:,2:end]))
    # Make suer we're ordered chronologically
    decomp_df = sort!(decomp_df, :year)

    decomp_vars = ["LE_mod", "H_mod", "h_mod", "Lstar_mod", "Lmed_mod",
        "LE_b", "LE_B", "LE_c", "LE_C", "LE_d",
        "H_b", "H_B", "H_c", "H_C", "H_d",
        "h_b", "h_B", "h_c", "h_C", "h_d",
        "Lstar_b", "Lstar_B", "Lstar_c", "Lstar_C", "Lstar_d",
        "Lmed_b", "Lmed_B", "Lmed_c", "Lmed_C", "Lmed_d",
        "Δb", "ΔB", "Δc", "ΔC", "Δd",
        "ΔLE_mod", "ΔH_mod", "Δh_mod", "ΔLstar_mod", "ΔLmed_mod"]
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
        decomp_df.h_mod[ii] = h(params, eval_age, spec = spec)
        decomp_df.Lstar_mod[ii] = lifespan(params, Sstar = 0.001, spec = spec)
        decomp_df.Lmed_mod[ii] = Lmed(params, eval_age, spec = spec)
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
        # Compute h gradient for each parameter
        decomp_df.h_B[ii] = hgrad(params, eval_age, :B, spec = spec)
        decomp_df.h_b[ii] = hgrad(params, eval_age, :b, spec = spec)
        decomp_df.h_C[ii] = hgrad(params, eval_age, :C, spec = spec)
        decomp_df.h_c[ii] = hgrad(params, eval_age, :c, spec = spec)
        decomp_df.h_d[ii] = hgrad(params, eval_age, :d, spec = spec)
        # Compute Lstar gradient for each parameter
        decomp_df.Lstar_B[ii] = lifespangrad(params, 0.001, :B, spec = spec)
        decomp_df.Lstar_b[ii] = lifespangrad(params, 0.001, :b, spec = spec)
        decomp_df.Lstar_C[ii] = lifespangrad(params, 0.001, :C, spec = spec)
        decomp_df.Lstar_c[ii] = lifespangrad(params, 0.001, :c, spec = spec)
        decomp_df.Lstar_d[ii] = lifespangrad(params, 0.001, :d, spec = spec)
        # Compute Lmed gradient for each parameter
        decomp_df.Lmed_B[ii] = Lmedgrad(params, eval_age, :B, spec = spec)
        decomp_df.Lmed_b[ii] = Lmedgrad(params, eval_age, :b, spec = spec)
        decomp_df.Lmed_C[ii] = Lmedgrad(params, eval_age, :C, spec = spec)
        decomp_df.Lmed_c[ii] = Lmedgrad(params, eval_age, :c, spec = spec)
        decomp_df.Lmed_d[ii] = Lmedgrad(params, eval_age, :d, spec = spec)

        if ii > 1
            decomp_df.ΔB[ii] = decomp_df.B[ii] - decomp_df.B[ii-1]
            decomp_df.Δb[ii] = decomp_df.b[ii] - decomp_df.b[ii-1]
            decomp_df.ΔC[ii] = decomp_df.C[ii] - decomp_df.C[ii-1]
            decomp_df.Δc[ii] = decomp_df.c[ii] - decomp_df.c[ii-1]
            decomp_df.Δd[ii] = decomp_df.d[ii] - decomp_df.d[ii-1]
            decomp_df.ΔLE_mod[ii] = decomp_df.LE_mod[ii] - decomp_df.LE_mod[ii-1]
            decomp_df.ΔH_mod[ii] = decomp_df.H_mod[ii] - decomp_df.H_mod[ii-1]
            decomp_df.Δh_mod[ii] = decomp_df.h_mod[ii] - decomp_df.h_mod[ii-1]
            decomp_df.ΔLstar_mod[ii] = decomp_df.Lstar_mod[ii] - decomp_df.Lstar_mod[ii-1]
            decomp_df.ΔLmed_mod[ii] = decomp_df.Lmed_mod[ii] - decomp_df.Lmed_mod[ii-1]
        end

    end

    insertcols!(decomp_df, 2, :forecast => repeat([0], nrow(decomp_df)) )

    if forecasts
        frc_years = unique(parests_df.year[parests_df.forecast.==1])
        decomp_df.forecast[in.(decomp_df.year, [frc_years])] .= 1
    end


    return decomp_df

end




"""
A function that extracts a vector of SilerParam structures from rows of a DataFrame
    of samples from a posterior distribution
"""
function particles2params(df_particle, Tptt; log_pars = true, model_vers = :indep)

    if  model_vers == :cov
        B = Symbol.("lpars[".*string.(Tptt).*"][1]")
        b = Symbol.("lpars[".*string.(Tptt).*"][2]")
        C = Symbol.("lpars[".*string.(Tptt).*"][3]")
        c = Symbol.("lpars[".*string.(Tptt).*"][4]")
        d = Symbol.("lpars[".*string.(Tptt).*"][5]")
    else
        B = Symbol("lB["*string(Tptt)*"]")
        b = Symbol("lb["*string(Tptt)*"]")
        C = Symbol("lC["*string(Tptt)*"]")
        c = Symbol("lc["*string(Tptt)*"]")
        d = Symbol("ld["*string(Tptt)*"]")
    end


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



"""
 Compute the model implied LE, H and h for the in-sample periods
"""
function compute_LE_post(df_post, years, nahead; spec = :Bergeron, model_vers = :indep)
    # Add extra columns for model impled LE and H
    LEs = Symbol.("LE[".*string.(1:length(years)+nahead).*"]")
    Hs = Symbol.("H[".*string.(1:length(years)+nahead).*"]")
    hs = Symbol.("h[".*string.(1:length(years)+nahead).*"]")
    Lstars_99p9 = Symbol.("Lstar_99p9[".*string.(1:length(years)+nahead).*"]")
    Lstars_99 = Symbol.("Lstar_99[".*string.(1:length(years)+nahead).*"]")
    Lstars_95 = Symbol.("Lstar_95[".*string.(1:length(years)+nahead).*"]")
    Lstars_90 = Symbol.("Lstar_90[".*string.(1:length(years)+nahead).*"]")
    Lmeds = Symbol.("Lmed[".*string.(1:length(years)+nahead).*"]")
    impl_vars = vcat(LEs, Hs, hs, Lstars_99p9, Lstars_99, Lstars_95, Lstars_90)
    df_post = hcat(df_post,DataFrame(NaN.*zeros(nrow(df_post), length(impl_vars)), impl_vars))
    # Loop through each period
    prog = Progress(length(years), desc = "Calculating model implied LE and H: ")
    for ii in 1:length(years)
        params = particles2params(df_post, ii, log_pars = true, model_vers = model_vers)
        df_post[:, LEs[ii]] = LE.(params, [0.0], spec = spec)
        df_post[:, Hs[ii]] = H.(params, [0.0], spec = spec)
        df_post[:, hs[ii]] = h.(params, [0.0], spec = spec)
        df_post[:, Lstars_99p9[ii]] = lifespan.(params, Sstar = 0.001, spec = spec)
        df_post[:, Lstars_99[ii]] = lifespan.(params, Sstar = 0.01, spec = spec)
        df_post[:, Lstars_95[ii]] = lifespan.(params, Sstar = 0.05, spec = spec)
        df_post[:, Lstars_90[ii]] = lifespan.(params, Sstar = 0.1, spec = spec)
        df_post[:, Lmeds[ii]] = Lmed.(params, [0.0], spec = spec)
        next!(prog)
    end

    return df_post

end



"""
Function to compute predictive distribution nahead periods ahead with ndraws draws for
    the path of future shocks
"""
function compute_forecasts(df_post, nahead, ndraws, years; spec = :Bergeron, model_vers = :indep)

    # Siler parameters
    if model_vers == :cov
        Bs = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][1]")
        bs = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][2]")
        Cs = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][3]")
        cs = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][4]")
        ds = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][5]")
        σs = Symbol.("lpars[".*string.(1:length(years)+nahead).*"][6]")
        # Future parameters
        B_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][1]")
        b_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][2]")
        C_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][3]")
        c_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][4]")
        d_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][5]")
        σ_f = Symbol.("lpars[".*string.(length(years)+1:length(years)+nahead).*"][6]")
    else
        Bs = Symbol.("lB[".*string.(1:length(years)+nahead).*"]")
        bs = Symbol.("lb[".*string.(1:length(years)+nahead).*"]")
        Cs = Symbol.("lC[".*string.(1:length(years)+nahead).*"]")
        cs = Symbol.("lc[".*string.(1:length(years)+nahead).*"]")
        ds = Symbol.("ld[".*string.(1:length(years)+nahead).*"]")
        σs = Symbol.("lσ[".*string.(1:length(years)+nahead).*"]")
        # Future parameters
        B_f = Symbol.("lB[".*string.(length(years)+1:length(years)+nahead).*"]")
        b_f = Symbol.("lb[".*string.(length(years)+1:length(years)+nahead).*"]")
        C_f = Symbol.("lC[".*string.(length(years)+1:length(years)+nahead).*"]")
        c_f = Symbol.("lc[".*string.(length(years)+1:length(years)+nahead).*"]")
        d_f = Symbol.("ld[".*string.(length(years)+1:length(years)+nahead).*"]")
        σ_f = Symbol.("lσ[".*string.(length(years)+1:length(years)+nahead).*"]")
    end
    # Drift terms (remeber lB[tt] = α_B[tt-1] + lB[tt-1] = shock)
    α_Bs = Symbol.("α_B[".*string.(1:length(years)-1+nahead).*"]")
    α_bs = Symbol.("α_b[".*string.(1:length(years)-1+nahead).*"]")
    α_Cs = Symbol.("α_C[".*string.(1:length(years)-1+nahead).*"]")
    α_cs = Symbol.("α_c[".*string.(1:length(years)-1+nahead).*"]")
    α_ds = Symbol.("α_d[".*string.(1:length(years)-1+nahead).*"]")
    α_σs = Symbol.("α_σ[".*string.(1:length(years)-1+nahead).*"]")
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
    hs = Symbol.("h[".*string.(1:length(years)+nahead).*"]")
    Lstars_99p9 = Symbol.("Lstar_99p9[".*string.(1:length(years)+nahead).*"]")
    Lstars_99 = Symbol.("Lstar_99[".*string.(1:length(years)+nahead).*"]")
    Lstars_95 = Symbol.("Lstar_95[".*string.(1:length(years)+nahead).*"]")
    Lstars_90 = Symbol.("Lstar_90[".*string.(1:length(years)+nahead).*"]")
    Lmeds = Symbol.("Lmed[".*string.(1:length(years)+nahead).*"]")

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
            if model_vers == :cov
                ρ_bB = df_particle[1,Symbol("ρ_bB")]
                ρ_cC = df_particle[1,Symbol("ρ_cC")]
                ρ_cd = df_particle[1,Symbol("ρ_cd")]
                ρ_Cd = df_particle[1,Symbol("ρ_Cd")]
            end

            # Forecast drift terms (remeber that we need to go one further back for these)
            df_particle[:,α_Bs[Tptt-1]] = df_particle[:,α_Bs[Tptt-2]] + rand(Normal(0.0, σ_ξB),ndraws)
            df_particle[:,α_bs[Tptt-1]] = df_particle[:,α_bs[Tptt-2]] + rand(Normal(0.0, σ_ξb),ndraws)
            df_particle[:,α_Cs[Tptt-1]] = df_particle[:,α_Cs[Tptt-2]] + rand(Normal(0.0, σ_ξC),ndraws)
            df_particle[:,α_cs[Tptt-1]] = df_particle[:,α_cs[Tptt-2]] + rand(Normal(0.0, σ_ξc),ndraws)
            df_particle[:,α_ds[Tptt-1]] = df_particle[:,α_ds[Tptt-2]] + rand(Normal(0.0, σ_ξd),ndraws)
            df_particle[:,α_σs[Tptt-1]] = df_particle[:,α_σs[Tptt-2]] + rand(Normal(0.0, σ_ξσ),ndraws)

            # Forecast the parameters
            if model_vers == :cov
                # Extract the covariance matrix
                σ_pars = [σ_ϵB, σ_ϵb, σ_ϵC, σ_ϵc, σ_ϵd, σ_ϵσ]
                Σ_ϵ = Matrix(Diagonal(σ_pars))
                Σ_ϵ[1,2] = ρ_bB*sqrt(σ_ϵB)*sqrt(σ_ϵb)
                Σ_ϵ[2,1] = ρ_bB*sqrt(σ_ϵB)*sqrt(σ_ϵb)
                Σ_ϵ[3,4] = ρ_cC*sqrt(σ_ϵC)*sqrt(σ_ϵc)
                Σ_ϵ[4,3] = ρ_cC*sqrt(σ_ϵC)*sqrt(σ_ϵc)
                Σ_ϵ[3,5] = ρ_Cd*sqrt(σ_ϵC)*sqrt(σ_ϵd)
                Σ_ϵ[5,3] = ρ_Cd*sqrt(σ_ϵC)*sqrt(σ_ϵd)
                Σ_ϵ[4,5] = ρ_cd*sqrt(σ_ϵc)*sqrt(σ_ϵd)
                Σ_ϵ[5,4] = ρ_cd*sqrt(σ_ϵc)*sqrt(σ_ϵd)
                # Draw future shocks
                ϵ_draws = rand(MvNormal(zeros(6), Σ_ϵ),ndraws)
                # Update the parameters
                df_particle[:,Bs[Tptt]] = df_particle[:,α_Bs[Tptt-1]] + df_particle[:,Bs[Tptt-1]] + ϵ_draws[1,:]
                df_particle[:,bs[Tptt]] = df_particle[:,α_bs[Tptt-1]] + df_particle[:,bs[Tptt-1]] + ϵ_draws[2,:]
                df_particle[:,Cs[Tptt]] = df_particle[:,α_Cs[Tptt-1]] + df_particle[:,Cs[Tptt-1]] + ϵ_draws[3,:]
                df_particle[:,cs[Tptt]] = df_particle[:,α_cs[Tptt-1]] + df_particle[:,cs[Tptt-1]] + ϵ_draws[4,:]
                df_particle[:,ds[Tptt]] = df_particle[:,α_ds[Tptt-1]] + df_particle[:,ds[Tptt-1]] + ϵ_draws[5,:]
                df_particle[:,σs[Tptt]] = df_particle[:,α_σs[Tptt-1]] + df_particle[:,σs[Tptt-1]] + ϵ_draws[6,:]
            else
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
            end

            # Compute the model implied LE and H overtime
            params = particles2params(df_particle, Tptt, log_pars = true, model_vers = model_vers)
            df_particle[:,LEs[Tptt]] = LE.(params, [0.0], spec = spec)
            df_particle[:,Hs[Tptt]] = H.(params, [0.0], spec = spec)
            df_particle[:,hs[Tptt]] = h.(params, [0.0], spec = spec)
            df_particle[:, Lstars_99p9[Tptt]] = lifespan.(params, Sstar = 0.001, spec = spec)
            df_particle[:, Lstars_99[Tptt]] = lifespan.(params, Sstar = 0.01, spec = spec)
            df_particle[:, Lstars_95[Tptt]] = lifespan.(params, Sstar = 0.05, spec = spec)
            df_particle[:, Lstars_90[Tptt]] = lifespan.(params, Sstar = 0.1, spec = spec)
            #df_particle[:,Lstars[Tptt]] = lifespan.(params, Sstar = 0.001, spec = spec)
            df_particle[:,Lmeds[Tptt]] = Lmed.(params, [0.0], spec = spec)

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
        log_pars = true, spec = :Bergeron, model_vers = :i2drift)


    # Work on deep copy version as we might transform to account for logs
    df_in = deepcopy(df_pred)
    if (model_vers != :i2drift) & (model_vers != :cov)
        throw("Only supported for i2drift and cov")
    end
    all_years = vcat(past_years, fut_years)
    # Siler parameter names
    if log_pars
        if model_vers == :cov
            Bs = Symbol.("lpars[".*string.(1:length(all_years)).*"][1]")
            bs = Symbol.("lpars[".*string.(1:length(all_years)).*"][2]")
            Cs = Symbol.("lpars[".*string.(1:length(all_years)).*"][3]")
            cs = Symbol.("lpars[".*string.(1:length(all_years)).*"][4]")
            ds = Symbol.("lpars[".*string.(1:length(all_years)).*"][5]")
            σs = Symbol.("lpars[".*string.(1:length(all_years)).*"][6]")
        else
            Bs = Symbol.("lB[".*string.(1:length(all_years)).*"]")
            bs = Symbol.("lb[".*string.(1:length(all_years)).*"]")
            Cs = Symbol.("lC[".*string.(1:length(all_years)).*"]")
            cs = Symbol.("lc[".*string.(1:length(all_years)).*"]")
            ds = Symbol.("ld[".*string.(1:length(all_years)).*"]")
            σs = Symbol.("lσ[".*string.(1:length(all_years)).*"]")
        end
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
        df_in[:, Bs] = Matrix(df_in[:, bs]).*Matrix(df_in[:, Bs]) .- log.(Matrix(df_in[:, bs]))
        df_in[:, Cs] = Matrix(df_in[:, cs]).*Matrix(df_in[:, Cs]) .- log.(Matrix(df_in[:, cs]))
    elseif spec == :Scott
        df_in[:, Bs] = Matrix(df_in[:, Bs]) .- log.(Matrix(df_in[:, bs]))./Matrix(df_in[:, bs])
        df_in[:, Cs] = Matrix(df_in[:, Cs]) .- log.(Matrix(df_in[:, cs]))./Matrix(df_in[:, cs])
    elseif spec == :Bergeron
        df_in[:, Bs] = df_in[:, Bs]
        df_in[:, Cs] = df_in[:, Cs]
    elseif spec == :Standard
        throw("Haven't worked out the new conversion for Standard yet")
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

    elseif (model_vers == :i2drift) | (model_vers == :cov)

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

        if model_vers == :cov
            ρ_s = [:ρ_bB, :ρ_cC, :ρ_cd, :ρ_Cd]
            ρ_ests = summarise_forecasts(df_in[:,ρ_s], Int.(1:length(ρ_s)), parname = :ρ)
            ρ_ests.parameter = [:ρ_bB, :ρ_cC, :ρ_cd, :ρ_Cd]
            ρ_ests.year .= 0
            par_ests = vcat(par_ests, ρ_ests)
        end

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
    hs = Symbol.("h[".*string.(1:length(all_years)).*"]")
    Lstars_99p9 = Symbol.("Lstar_99p9[".*string.(1:length(all_years)).*"]")
    Lstars_99 = Symbol.("Lstar_99[".*string.(1:length(all_years)).*"]")
    Lstars_95 = Symbol.("Lstar_95[".*string.(1:length(all_years)).*"]")
    Lstars_90 = Symbol.("Lstar_90[".*string.(1:length(all_years)).*"]")
    #Lstars = Symbol.("Lstar[".*string.(1:length(all_years)).*"]")
    Lmeds = Symbol.("Lmed[".*string.(1:length(all_years)).*"]")


    LE_ests = summarise_forecasts(df_in[:,LEs], all_years, parname = :LE)
    H_ests = summarise_forecasts(df_in[:,Hs], all_years, parname = :H)
    h_ests = summarise_forecasts(df_in[:,hs], all_years, parname = :h)
    Lstar_99p9_ests = summarise_forecasts(df_in[:,Lstars_99p9], all_years, parname = :Lstar_99p9)
    Lstar_99_ests = summarise_forecasts(df_in[:,Lstars_99], all_years, parname = :Lstar_99)
    Lstar_95_ests = summarise_forecasts(df_in[:,Lstars_95], all_years, parname = :Lstar_95)
    Lstar_90_ests = summarise_forecasts(df_in[:,Lstars_90], all_years, parname = :Lstar_90)
    #Lstar_ests = summarise_forecasts(df_in[:,Lstars], all_years, parname = :Lstar)
    Lmed_ests = summarise_forecasts(df_in[:,Lmeds], all_years, parname = :Lmed)


    par_ests = vcat(par_ests, LE_ests, H_ests, h_ests, Lstar_99p9_ests, Lstar_99_ests,
        Lstar_95_ests, Lstar_90_ests, Lmed_ests)


    insertcols!(par_ests, 2, :forecast => repeat([0], nrow(par_ests)) )
    par_ests.forecast[in.(par_ests.year, [vcat(fut_years, past_years[end])])] .= 1

    return par_ests

end









"""
Function to compute gradients across the age distribution from a decomp_df
"""
function compute_LEgrad_df(decomp_pred; ages = 0:140)
    # Years
    all_years = decomp_pred.year
    # Initialise dataframe for results
    LEgrad_df = DataFrame(age = repeat(ages, length(all_years)), year = repeat(all_years, inner = length(ages)),
        forecast = repeat(decomp_pred.forecast, inner = length(ages)))
    grad_vars = [:LE, :LE_Bs, :LE_bs, :LE_Cs, :LE_cs, :LE_ds, :LE_cC,
        :Lstar, :Lstar_Bs, :Lstar_bs, :Lstar_Cs, :Lstar_cs, :Lstar_ds, :Lstar_cC,
        :Lmed, :Lmed_Bs, :Lmed_bs, :Lmed_Cs, :Lmed_cs, :Lmed_ds, :Lmed_cC,
        :h, :h_Bs, :h_bs, :h_Cs, :h_cs, :h_ds, :h_cC,
        :H, :H_Bs, :H_bs, :H_Cs, :H_cs, :H_ds, :H_cC,
        :mortality, :survival]
    LEgrad_df = hcat(LEgrad_df,DataFrame(NaN.*zeros(nrow(LEgrad_df), length(grad_vars)), grad_vars))
    # Loop over each year
    prog = Progress(length(all_years), desc = "Computing gradient df: ")
    for yy in all_years
        # Get parameters for this year
        obs = Int.(1:nrow(decomp_pred))[decomp_pred.year .== yy][1]
        row = NamedTuple(decomp_pred[obs,:])
        params = SilerParam(b = row.b, B = row.B, c = row.c, C = row.C, d = row.d)

        # Fill grad dataframe
        out_obs = Int.(1:nrow(LEgrad_df))[LEgrad_df.year .== yy]
        # LE
        LEgrad_df.LE[out_obs] = LE.([params], ages, spec = :Bergeron)
        LEgrad_df.LE_Bs[out_obs] = LEgrad.([params], ages, [:B], spec = :Bergeron)
        LEgrad_df.LE_bs[out_obs] = LEgrad.([params], ages, [:b], spec = :Bergeron)
        LEgrad_df.LE_Cs[out_obs] = LEgrad.([params], ages, [:C], spec = :Bergeron)
        LEgrad_df.LE_cs[out_obs] = LEgrad.([params], ages, [:c], spec = :Bergeron)
        LEgrad_df.LE_ds[out_obs] = LEgrad.([params], ages, [:d], spec = :Bergeron)
        LEgrad_df.LE_cC[out_obs] = LEcross.([params], ages, [:c], [:C], spec = :Bergeron)
        # Lstar
        LEgrad_df.Lstar[out_obs] .= lifespan.([params], Sstar = 0.001, spec = :Bergeron)
        LEgrad_df.Lstar_Bs[out_obs] .= lifespangrad.([params], [0.001], [:B], spec = :Bergeron)
        LEgrad_df.Lstar_bs[out_obs] .= lifespangrad.([params], [0.001], [:b], spec = :Bergeron)
        LEgrad_df.Lstar_Cs[out_obs] .= lifespangrad.([params], [0.001], [:C], spec = :Bergeron)
        LEgrad_df.Lstar_cs[out_obs] .= lifespangrad.([params], [0.001], [:c], spec = :Bergeron)
        LEgrad_df.Lstar_ds[out_obs] .= lifespangrad.([params], [0.001], [:d], spec = :Bergeron)
        # Lmed
        LEgrad_df.Lmed[out_obs] = Lmed.([params], ages, spec = :Bergeron)
        LEgrad_df.Lmed_Bs[out_obs] = Lmedgrad.([params], ages, [:B], spec = :Bergeron)
        LEgrad_df.Lmed_bs[out_obs] = Lmedgrad.([params], ages, [:b], spec = :Bergeron)
        LEgrad_df.Lmed_Cs[out_obs] = Lmedgrad.([params], ages, [:C], spec = :Bergeron)
        LEgrad_df.Lmed_cs[out_obs] = Lmedgrad.([params], ages, [:c], spec = :Bergeron)
        LEgrad_df.Lmed_ds[out_obs] = Lmedgrad.([params], ages, [:d], spec = :Bergeron)
        # H
        LEgrad_df.H[out_obs] = H.([params], ages, spec = :Bergeron)
        LEgrad_df.H_Bs[out_obs] = Hgrad.([params], ages, [:B], spec = :Bergeron)
        LEgrad_df.H_bs[out_obs] = Hgrad.([params], ages, [:b], spec = :Bergeron)
        LEgrad_df.H_Cs[out_obs] = Hgrad.([params], ages, [:C], spec = :Bergeron)
        LEgrad_df.H_cs[out_obs] = Hgrad.([params], ages, [:c], spec = :Bergeron)
        LEgrad_df.H_ds[out_obs] = Hgrad.([params], ages, [:d], spec = :Bergeron)
        # H
        LEgrad_df.h[out_obs] = h.([params], ages, spec = :Bergeron)
        LEgrad_df.h_Bs[out_obs] = hgrad.([params], ages, [:B], spec = :Bergeron)
        LEgrad_df.h_bs[out_obs] = hgrad.([params], ages, [:b], spec = :Bergeron)
        LEgrad_df.h_Cs[out_obs] = hgrad.([params], ages, [:C], spec = :Bergeron)
        LEgrad_df.h_cs[out_obs] = hgrad.([params], ages, [:c], spec = :Bergeron)
        LEgrad_df.h_ds[out_obs] = hgrad.([params], ages, [:d], spec = :Bergeron)
        LEgrad_df.h_cC[out_obs] = hcross.([params], ages, [:c], [:C], spec = :Bergeron)


        LEgrad_df.mortality[out_obs] = siler.([params], ages; spec = :Bergeron)
        LEgrad_df.survival[out_obs] = siler_S.([params], [0.0], ages; spec = :Bergeron)
        next!(prog)
    end

    return LEgrad_df

end
