"""
Helper functions to analyse Turing results
"""

"""
Function to extract summary statistics from time series variables
"""
function summarise_stats(df_post, years_selected; log_pars = false, parname = :unknown)

    if log_pars
        df_post[:, 3:end] = exp.(df_post[:, 3:end])
    end
    rename!(df_post, vcat(names(df_post)[1:2], string.(years_selected)))
    df_ests = describe(df_post[:, 3:end])
    rename!(df_ests, Dict(:variable => :year))
    df_ests.year = parse.(Int64,string.(df_ests.year))
    insertcols!(df_ests, 1, :parameter => repeat([parname], nrow(df_ests)) )

    df_ests[:,:pc975] .= 0.0
    df_ests[:,:pc025] .= 0.0
    df_ests[:,:pc85] .= 0.0
    df_ests[:,:pc15] .= 0.0
    for yy in 3:ncol(df_post)
        df_ests.pc975[(yy-2)] = percentile(df_post[:,yy], 97.5)
        df_ests.pc025[(yy-2)] = percentile(df_post[:,yy], 2.5)
        df_ests.pc85[(yy-2)] = percentile(df_post[:,yy], 85)
        df_ests.pc15[(yy-2)] = percentile(df_post[:,yy], 15)
    end

    return(df_ests)
end


"""
Function that creates summary dataframes for each Siler parameter
"""
function extract_variables(chain_in, periods, years_selected; log_pars = false)
    # Names
    if log_pars
        Bs = vcat([:iteration, :chain], Symbol.("lB[".*string.(1:length(periods)).*"]"))
        bs = vcat([:iteration, :chain], Symbol.("lb[".*string.(1:length(periods)).*"]"))
        Cs = vcat([:iteration, :chain], Symbol.("lC[".*string.(1:length(periods)).*"]"))
        cs = vcat([:iteration, :chain], Symbol.("lc[".*string.(1:length(periods)).*"]"))
        ds = vcat([:iteration, :chain], Symbol.("ld[".*string.(1:length(periods)).*"]"))
        σs = vcat([:iteration, :chain], Symbol.("lσ[".*string.(1:length(periods)).*"]"))
    else
        Bs = vcat([:iteration, :chain], Symbol.("B[".*string.(1:length(periods)).*"]"))
        bs = vcat([:iteration, :chain], Symbol.("b[".*string.(1:length(periods)).*"]"))
        Cs = vcat([:iteration, :chain], Symbol.("C[".*string.(1:length(periods)).*"]"))
        cs = vcat([:iteration, :chain], Symbol.("c[".*string.(1:length(periods)).*"]"))
        ds = vcat([:iteration, :chain], Symbol.("d[".*string.(1:length(periods)).*"]"))
        σs = vcat([:iteration, :chain], Symbol.("σ[".*string.(1:length(periods)).*"]"))
    end
    # Extract results
    df_indep = DataFrame(chain_in)
    B_ests = summarise_stats(df_indep[:,Bs], years_selected, log_pars = log_pars, parname = :B)
    b_ests = summarise_stats(df_indep[:,bs], years_selected, log_pars = log_pars, parname = :b)
    C_ests = summarise_stats(df_indep[:,Cs], years_selected, log_pars = log_pars, parname = :C)
    c_ests = summarise_stats(df_indep[:,cs], years_selected, log_pars = log_pars, parname = :c)
    d_ests = summarise_stats(df_indep[:,ds], years_selected, log_pars = log_pars, parname = :d)
    σ_ests = summarise_stats(df_indep[:,σs], years_selected, log_pars = log_pars, parname = :σ)


    par_ests = vcat(B_ests, b_ests, C_ests, c_ests, d_ests, σ_ests)

    return par_ests

end


"""
Plot each Siler parameter over time
"""
function plot_siler_params(par_ests::DataFrame)

    plt = plot(layout = (2,3), xrotation = 45.0, margin=3Plots.mm, size = (800,400))

    B_ests = par_ests[par_ests.parameter .== :B,:]
    plot!(B_ests.year, B_ests.median, title = L"B_{t}", label = false, color = 1, subplot = 1)
    plot!(B_ests.year, B_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(B_ests.year, B_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(B_ests.year, B_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
    plot!(B_ests.year, B_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)

    b_ests = par_ests[par_ests.parameter .== :b,:]
    plot!(b_ests.year, b_ests.median, title = L"b_{t}", label = false, color = 1, subplot = 2)
    plot!(b_ests.year, b_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(b_ests.year, b_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(b_ests.year, b_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
    plot!(b_ests.year, b_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)

    C_ests = par_ests[par_ests.parameter .== :C,:]
    plot!(C_ests.year, C_ests.median, title = L"C_{t}", label = false, color = 1, subplot = 4)
    plot!(C_ests.year, C_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(C_ests.year, C_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(C_ests.year, C_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 4)
    plot!(C_ests.year, C_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 4)

    c_ests = par_ests[par_ests.parameter .== :c,:]
    plot!(c_ests.year, c_ests.median, title = L"c_{t}", label = false, color = 1, subplot = 5)
    plot!(c_ests.year, c_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(c_ests.year, c_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(c_ests.year, c_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 5)
    plot!(c_ests.year, c_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 5)

    d_ests = par_ests[par_ests.parameter .== :d,:]
    plot!(d_ests.year, d_ests.median, title = L"d_{t}", label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 3)

    σ_ests = par_ests[par_ests.parameter .== :σ,:]
    plot!(σ_ests.year, σ_ests.median, title = L"\sigma_{t}", label = false, color = 1, subplot = 6)#, ylim = (0.0, 0.1))
    plot!(σ_ests.year, σ_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 6)

    return plt

end
