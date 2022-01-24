"""
Helper functions to analyse Turing results
"""

"""
Function to extract summary statistics from time series variables
"""
function summarise_stats(df_post, years_selected)
    rename!(df_post, vcat(names(df_post)[1:2], string.(years_selected)))
    df_ests = describe(df_post[:, 3:end])
    df_ests.variable = parse.(Int64,string.(df_ests.variable))
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
function extract_variables(chain_in, periods, years_selected)
    # Names
    Bs = vcat([:iteration, :chain], Symbol.("B[".*string.(1:length(periods)).*"]"))
    bs = vcat([:iteration, :chain], Symbol.("b[".*string.(1:length(periods)).*"]"))
    Cs = vcat([:iteration, :chain], Symbol.("C[".*string.(1:length(periods)).*"]"))
    cs = vcat([:iteration, :chain], Symbol.("c[".*string.(1:length(periods)).*"]"))
    ds = vcat([:iteration, :chain], Symbol.("d[".*string.(1:length(periods)).*"]"))
    σs = vcat([:iteration, :chain], Symbol.("σ[".*string.(1:length(periods)).*"]"))
    # Extract results
    df_indep = DataFrame(chain_in)
    B_ests = summarise_stats(df_indep[:,Bs], years_selected)
    b_ests = summarise_stats(df_indep[:,bs], years_selected)
    C_ests = summarise_stats(df_indep[:,Cs], years_selected)
    c_ests = summarise_stats(df_indep[:,cs], years_selected)
    d_ests = summarise_stats(df_indep[:,ds], years_selected)
    σ_ests = summarise_stats(df_indep[:,σs], years_selected)


    par_ests = (B=B_ests, b=b_ests, C=C_ests, c=c_ests, d=d_ests, σ=σ_ests)
    return par_ests

end


"""
Plot each Siler parameter over time
"""
function plot_siler_params(par_ests)

    plt = plot(layout = (2,3), xrotation = 45.0, margin=3Plots.mm, size = (800,400))
    plot!(par_ests.B.variable, par_ests.B.median, title = L"B_{t}", label = false, color = 1, subplot = 1)
    plot!(par_ests.B.variable, par_ests.B.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(par_ests.B.variable, par_ests.B.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(par_ests.B.variable, par_ests.B.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
    plot!(par_ests.B.variable, par_ests.B.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)

    plot!(par_ests.b.variable, par_ests.b.median, title = L"b_{t}", label = false, color = 1, subplot = 2)
    plot!(par_ests.b.variable, par_ests.b.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(par_ests.b.variable, par_ests.b.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(par_ests.b.variable, par_ests.b.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
    plot!(par_ests.b.variable, par_ests.b.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)

    plot!(par_ests.C.variable, par_ests.C.median, title = L"C_{t}", label = false, color = 1, subplot = 4)
    plot!(par_ests.C.variable, par_ests.C.pc975, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(par_ests.C.variable, par_ests.C.pc025, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(par_ests.C.variable, par_ests.C.pc85, linestyle = :dash, label = false, color = 1, subplot = 4)
    plot!(par_ests.C.variable, par_ests.C.pc15, linestyle = :dash, label = false, color = 1, subplot = 4)

    plot!(par_ests.c.variable, par_ests.c.median, title = L"c_{t}", label = false, color = 1, subplot = 5)
    plot!(par_ests.c.variable, par_ests.c.pc975, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(par_ests.c.variable, par_ests.c.pc025, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(par_ests.c.variable, par_ests.c.pc85, linestyle = :dash, label = false, color = 1, subplot = 5)
    plot!(par_ests.c.variable, par_ests.c.pc15, linestyle = :dash, label = false, color = 1, subplot = 5)

    plot!(par_ests.d.variable, par_ests.d.median, title = L"d_{t}", label = false, color = 1, subplot = 3)
    plot!(par_ests.d.variable, par_ests.d.pc975, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(par_ests.d.variable, par_ests.d.pc025, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(par_ests.d.variable, par_ests.d.pc85, linestyle = :dash, label = false, color = 1, subplot = 3)
    plot!(par_ests.d.variable, par_ests.d.pc15, linestyle = :dash, label = false, color = 1, subplot = 3)

    plot!(par_ests.σ.variable, par_ests.σ.median, title = L"\sigma_{t}", label = false, color = 1, subplot = 6)#, ylim = (0.0, 0.1))
    plot!(par_ests.σ.variable, par_ests.σ.pc975, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(par_ests.σ.variable, par_ests.σ.pc025, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(par_ests.σ.variable, par_ests.σ.pc85, linestyle = :dash, label = false, color = 1, subplot = 6)
    plot!(par_ests.σ.variable, par_ests.σ.pc15, linestyle = :dash, label = false, color = 1, subplot = 6)

    return plt

end
