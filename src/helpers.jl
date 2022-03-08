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
Basic Siler function to plot from paramters
"""
function siler(B,b,C,c,d, ages)

    μ = exp.(- b.* (ages .+ B)) .+ exp.(c .* (ages.- C)) .+ d
    lmort = log.(μ)
    mort = exp.(lmort)

    return mort
end



"""
Function to extract summary statistics from time series variables
"""
function summarise_stats(df_post, stats_post, years; log_pars = false, parname = :unknown)

    if log_pars
        df_post[:, 3:end] = exp.(df_post[:, 3:end])
    end
    rename!(df_post, vcat(names(df_post)[1:2], string.(years)))
    df_ests = describe(df_post[:, 3:end])
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
    for yy in 3:ncol(df_post)
        df_ests.std[(yy-2)] = std(df_post[:,yy])
        df_ests.pc975[(yy-2)] = percentile(df_post[:,yy], 97.5)
        df_ests.pc025[(yy-2)] = percentile(df_post[:,yy], 2.5)
        df_ests.pc85[(yy-2)] = percentile(df_post[:,yy], 85)
        df_ests.pc15[(yy-2)] = percentile(df_post[:,yy], 15)
        df_ests.pc75[(yy-2)] = percentile(df_post[:,yy], 75)
        df_ests.pc25[(yy-2)] = percentile(df_post[:,yy], 25)
    end
    # Add some convergence statistics
    df_out = hcat(df_ests, stats_post[:,[:mcse, :ess, :rhat, :ess_per_sec]])

    return(df_out)
end


"""
Function that creates summary dataframes for each Siler parameter
"""
function extract_variables(chain_in, years; log_pars = false, σ_pars = true, ext = false)
    # Names
    if log_pars
        if ext
            Bs = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][1]"))
            bs = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][2]"))
            Cs = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][3]"))
            cs = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][4]"))
            ds = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][5]"))
            σs = vcat([:iteration, :chain], Symbol.("pars[:,".*string.(1:length(years)).*"][6]"))
        else
            Bs = vcat([:iteration, :chain], Symbol.("lB[".*string.(1:length(years)).*"]"))
            bs = vcat([:iteration, :chain], Symbol.("lb[".*string.(1:length(years)).*"]"))
            Cs = vcat([:iteration, :chain], Symbol.("lC[".*string.(1:length(years)).*"]"))
            cs = vcat([:iteration, :chain], Symbol.("lc[".*string.(1:length(years)).*"]"))
            ds = vcat([:iteration, :chain], Symbol.("ld[".*string.(1:length(years)).*"]"))
            σs = vcat([:iteration, :chain], Symbol.("lσ[".*string.(1:length(years)).*"]"))
        end
    else
        Bs = vcat([:iteration, :chain], Symbol.("B[".*string.(1:length(years)).*"]"))
        bs = vcat([:iteration, :chain], Symbol.("b[".*string.(1:length(years)).*"]"))
        Cs = vcat([:iteration, :chain], Symbol.("C[".*string.(1:length(years)).*"]"))
        cs = vcat([:iteration, :chain], Symbol.("c[".*string.(1:length(years)).*"]"))
        ds = vcat([:iteration, :chain], Symbol.("d[".*string.(1:length(years)).*"]"))
        σs = vcat([:iteration, :chain], Symbol.("σ[".*string.(1:length(years)).*"]"))
    end
    # Extract results
    df_indep = DataFrame(chain_in)
    stats = DataFrame(summarystats(chain_in))
    B_ests = summarise_stats(df_indep[:,Bs], stats[in.(stats.parameters, [Bs]),:],
        years, log_pars = log_pars, parname = :B)
    b_ests = summarise_stats(df_indep[:,bs], stats[in.(stats.parameters, [bs]),:],
        years, log_pars = log_pars, parname = :b)
    C_ests = summarise_stats(df_indep[:,Cs], stats[in.(stats.parameters, [Cs]),:],
        years, log_pars = log_pars, parname = :C)
    c_ests = summarise_stats(df_indep[:,cs], stats[in.(stats.parameters, [cs]),:],
        years, log_pars = log_pars, parname = :c)
    d_ests = summarise_stats(df_indep[:,ds], stats[in.(stats.parameters, [ds]),:],
        years, log_pars = log_pars, parname = :d)
    σ_ests = summarise_stats(df_indep[:,σs], stats[in.(stats.parameters, [σs]),:],
        years, log_pars = log_pars, parname = :σ)

    par_ests = vcat(B_ests, b_ests, C_ests, c_ests, d_ests, σ_ests)

    # If we have a dynamic model with variance terms for the parameters, add these
    if σ_pars
        if ext
            # Parameter time series variance
            σ_par_names = vcat([:iteration, :chain], Symbol.("σ_par[".*string.(1:6).*"]"))
            σ_pars_ests = summarise_stats(df_indep[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_par)
            σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
            σ_pars_ests.year .= 0
            par_ests = vcat(par_ests, σ_pars_ests)
            # Parameter time series constant
            α_par_names = vcat([:iteration, :chain], Symbol.("α_pars[".*string.(1:6).*"]"))
            α_pars_ests = summarise_stats(df_indep[:,α_par_names], stats[in.(stats.parameters, [α_par_names]),:],
                Int.(1:6), log_pars = false, parname = :α_par)
            α_pars_ests.parameter = [:α_B, :α_b, :α_C, :α_c, :α_d, :α_σ]
            α_pars_ests.year .= 0
            par_ests = vcat(par_ests, α_pars_ests)
            # Parameter time series AR(1) coefficient
            β_par_names = vcat([:iteration, :chain], Symbol.("β_pars[".*string.(1:6).*"]"))
            β_pars_ests = summarise_stats(df_indep[:,β_par_names], stats[in.(stats.parameters, [β_par_names]),:],
                Int.(1:6), log_pars = false, parname = :β_par)
            β_pars_ests.parameter = [:β_B, :β_b, :β_C, :β_c, :β_d, :β_σ]
            β_pars_ests.year .= 0
            par_ests = vcat(par_ests, β_pars_ests)
            # Parameter time series linear trend
            τ_par_names = vcat([:iteration, :chain], Symbol.("τ_pars[".*string.(1:6).*"]"))
            τ_pars_ests = summarise_stats(df_indep[:,τ_par_names], stats[in.(stats.parameters, [τ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :τ_par)
            τ_pars_ests.parameter = [:τ_B, :τ_b, :τ_C, :τ_c, :τ_d, :τ_σ]
            τ_pars_ests.year .= 0
            par_ests = vcat(par_ests, τ_pars_ests)

        else
            σ_par_names = vcat([:iteration, :chain], Symbol.("σ_pars[".*string.(1:6).*"]"))
            σ_pars_ests = summarise_stats(df_indep[:,σ_par_names], stats[in.(stats.parameters, [σ_par_names]),:],
                Int.(1:6), log_pars = false, parname = :σ_par)
            σ_pars_ests.parameter = [:σ_B, :σ_b, :σ_C, :σ_c, :σ_d, :σ_σ]
            σ_pars_ests.year .= 0

            par_ests = vcat(par_ests, σ_pars_ests)
        end
    end

    return par_ests

end
