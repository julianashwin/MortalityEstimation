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

    rename!(df_sum, vcat(names(df_sum)[1:2], string.(years)))
    df_ests = describe(df_sum[:, 3:end])
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
    for yy in 3:ncol(df_sum)
        df_ests.std[(yy-2)] = std(df_sum[:,yy])
        df_ests.pc975[(yy-2)] = percentile(df_sum[:,yy], 97.5)
        df_ests.pc025[(yy-2)] = percentile(df_sum[:,yy], 2.5)
        df_ests.pc85[(yy-2)] = percentile(df_sum[:,yy], 85)
        df_ests.pc15[(yy-2)] = percentile(df_sum[:,yy], 15)
        df_ests.pc75[(yy-2)] = percentile(df_sum[:,yy], 75)
        df_ests.pc25[(yy-2)] = percentile(df_sum[:,yy], 25)
    end
    # Add some convergence statistics
    df_out = hcat(df_ests, stats_post[:,[:mcse, :ess, :rhat, :ess_per_sec]])

    return(df_out)
end


"""
Function that creates summary dataframes for each Siler parameter
    spec defines which Siler specification we want (Colchero, Scott, Bergeron or Standard)
    model_vers defines which dynamic model (indep, just_rw, i2drift, firstdiff)
"""
function extract_variables(chain_in, years; log_pars = false, spec = :Colchero,
        model_vers = :indep)

    # Siler parameter names
    if log_pars
        Bs = vcat(Symbol.("lB[".*string.(1:length(years)).*"]"))
        bs = vcat(Symbol.("lb[".*string.(1:length(years)).*"]"))
        Cs = vcat(Symbol.("lC[".*string.(1:length(years)).*"]"))
        cs = vcat(Symbol.("lc[".*string.(1:length(years)).*"]"))
        ds = vcat(Symbol.("ld[".*string.(1:length(years)).*"]"))
        σs = vcat(Symbol.("lσ[".*string.(1:length(years)).*"]"))
    else
        Bs = vcat(Symbol.("B[".*string.(1:length(years)).*"]"))
        bs = vcat(Symbol.("b[".*string.(1:length(years)).*"]"))
        Cs = vcat(Symbol.("C[".*string.(1:length(years)).*"]"))
        cs = vcat(Symbol.("c[".*string.(1:length(years)).*"]"))
        ds = vcat(Symbol.("d[".*string.(1:length(years)).*"]"))
        σs = vcat(Symbol.("σ[".*string.(1:length(years)).*"]"))
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
    if model_vers == :Colchero
        df_post[:, Bs] = df_post[:, Bs]
        df_post[:, Cs] = df_post[:, Cs]
    elseif model_vers == :Scott
        df_post[:, Bs] = Matrix(df_post[:, Bs])./Matrix(df_post[:, bs])
        df_post[:, Cs] = Matrix(df_post[:, Cs])./Matrix(df_post[:, cs])
    elseif model_vers == :Bergeron
        df_post[:, Bs] = exp.(.- Matrix(df_post[:, Bs]))
        df_post[:, Cs] = (Matrix(df_post[:, Cs]) .+ log.(Matrix(df_post[:, cs])))./Matrix(df_post[:, cs])
    elseif model_vers == :Standard
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

    # If we have a dynamic model with variance terms for the parameters, add these
    if σ_pars
        if firstdiff
            # Parameter time series variance
            σ_par_names = vcat([:iteration, :chain], Symbol.("σ_par[".*string.(1:6).*"]"))
            σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_par)
            σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
            σ_pars_ests.year .= 0
            par_ests = vcat(par_ests, σ_pars_ests)
            # Parameter time series constant
            α_par_names = vcat([:iteration, :chain], Symbol.("α_pars[".*string.(1:6).*"]"))
            α_pars_ests = summarise_stats(df_post[:,α_par_names], stats[in.(stats.parameters, [α_par_names]),:],
                Int.(1:6), log_pars = false, parname = :α_par)
            α_pars_ests.parameter = [:α_B, :α_b, :α_C, :α_c, :α_d, :α_σ]
            α_pars_ests.year .= 0
            par_ests = vcat(par_ests, α_pars_ests)
            # Parameter time series AR(1) coefficient
            β_par_names = vcat([:iteration, :chain], Symbol.("β_pars[".*string.(1:6).*"]"))
            β_pars_ests = summarise_stats(df_post[:,β_par_names], stats[in.(stats.parameters, [β_par_names]),:],
                Int.(1:6), log_pars = false, parname = :β_par)
            β_pars_ests.parameter = [:β_B, :β_b, :β_C, :β_c, :β_d, :β_σ]
            β_pars_ests.year .= 0
            par_ests = vcat(par_ests, β_pars_ests)
            # Parameter time series linear trend
            #τ_par_names = vcat([:iteration, :chain], Symbol.("τ_pars[".*string.(1:6).*"]"))
            #τ_pars_ests = summarise_stats(df_indep[:,τ_par_names], stats[in.(stats.parameters, [τ_par_names]),:],
            #    Int.(1:6), log_pars = false, parname = :τ_par)
            #τ_pars_ests.parameter = [:τ_B, :τ_b, :τ_C, :τ_c, :τ_d, :τ_σ]
            #τ_pars_ests.year .= 0
            #par_ests = vcat(par_ests, τ_pars_ests)
        elseif ext
            σ_par_names = vcat([:iteration, :chain], Symbol.("σ_pars[".*string.(1:6).*"]"))
            σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_par)
            σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
            σ_pars_ests.year .= 0

            par_ests = vcat(par_ests, σ_pars_ests)

            σ_α_par_names = vcat([:iteration, :chain], Symbol.("σ_α".*["B", "b", "C", "c", "d", "σ"]))
            σ_α_pars_ests = summarise_stats(df_post[:,(σ_α_par_names)], stats[in.(stats.parameters, [σ_α_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_α_par)
            σ_α_pars_ests.parameter = [:σ_αB, :σ_αb, :σ_αC, :σ_αc, :σ_αd, :σ_ασ]
            σ_α_pars_ests.year .= 0

            par_ests = vcat(par_ests, σ_α_pars_ests)

            α_Bs = vcat([:iteration, :chain], Symbol.("α_B[".*string.(1:length(years)-1).*"]"))
            α_bs = vcat([:iteration, :chain], Symbol.("α_b[".*string.(1:length(years)-1).*"]"))
            α_Cs = vcat([:iteration, :chain], Symbol.("α_C[".*string.(1:length(years)-1).*"]"))
            α_cs = vcat([:iteration, :chain], Symbol.("α_c[".*string.(1:length(years)-1).*"]"))
            α_ds = vcat([:iteration, :chain], Symbol.("α_d[".*string.(1:length(years)-1).*"]"))
            α_σs = vcat([:iteration, :chain], Symbol.("α_σ[".*string.(1:length(years)-1).*"]"))

            α_B_ests = summarise_stats(df_post[:,α_Bs], stats[in.(stats.parameters, [α_Bs]),:],
                years[2:end], log_pars = false, parname = :α_B)
            α_b_ests = summarise_stats(df_post[:,α_bs], stats[in.(stats.parameters, [α_bs]),:],
                years[2:end], log_pars = false, parname = :α_b)
            α_C_ests = summarise_stats(df_post[:,α_Cs], stats[in.(stats.parameters, [α_Cs]),:],
                years[2:end], log_pars = false, parname = :α_C)
            α_c_ests = summarise_stats(df_post[:,α_cs], stats[in.(stats.parameters, [α_cs]),:],
                years[2:end], log_pars = false, parname = :α_c)
            α_d_ests = summarise_stats(df_post[:,α_ds], stats[in.(stats.parameters, [α_ds]),:],
                years[2:end], log_pars = false, parname = :α_d)
            α_σ_ests = summarise_stats(df_post[:,α_σs], stats[in.(stats.parameters, [α_σs]),:],
                years[2:end], log_pars = false, parname = :α_σ)

            par_ests = vcat(par_ests, α_B_ests, α_b_ests, α_C_ests, α_c_ests, α_d_ests, α_σ_ests)

        else
            σ_par_names = vcat([:iteration, :chain], Symbol.("σ_pars[".*string.(1:6).*"]"))
            σ_pars_ests = summarise_stats(df_post[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_par)
            σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
            σ_pars_ests.year .= 0

            par_ests = vcat(par_ests, σ_pars_ests)
        end
    end

    return par_ests

end
