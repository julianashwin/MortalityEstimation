"""
Plot μ, S, LE and H given siler parameters
"""
function base_plt(illus_df, param; yscale = (0.,0.))
    plt = plot(legend = :top, size = (600,400), right_margin=12Plots.mm, left_margin = 8Plots.mm,
        xlabel = "Age", ylabel = "Rate", ylims = (0.0, 1.05), xlims = (-2,maximum(illus_df.age)))
    plt = plot!(illus_df.age, illus_df.μ_base, label = "Mortality", color = 1)
    plt = plot!(illus_df.age, illus_df.S_base, label = "Survival", color = 2)
    plt = plot!(illus_df.age, illus_df.H_base, label = "Inequality", color = 3)
    plt = plot!([0.0], [0.0], label = "LE", color = 4, legend = :top)
    if yscale == (0.,0.)
        plt = plot!(twinx(), illus_df.LE_base, label = "", color = 4, legend = :top,
        ylabel = "Remaining LE", xlims = (-2,maximum(illus_df.age)))
    else
        plt = plot!(twinx(), illus_df.LE_base, label = "", color = 4, legend = :top,
        ylabel = "Remaining LE", xlims = (-2,maximum(illus_df.age)), ylims = yscale)
    end

    return plt
end


"""
Plot changes to μ, S, LE and H driven by θ
"""
function param_change_plt(illus_df, param, θ; spec = :Colchero, LE_change = 10)
    LE_base = LE(param, 0.0, spec = spec)
    LE_new = LE_base + LE_change

    θ_new = find_zero(function(x) temp_param=deepcopy(param);
        setfield!(temp_param, θ, max(x, 1e-16));
        LE_new - LE(temp_param, 0.0, spec = spec) end,
        getfield(param, θ))
        param_new = deepcopy(param)
        setfield!(param_new, θ, θ_new)

    illus_df[:,Symbol("μ_"*string(θ)*"change")] = siler.([param_new], illus_df.age, spec = spec)
    illus_df[:,Symbol("S_"*string(θ)*"change")] = siler_S.([param_new], [0.0], illus_df.age, spec = spec)
    illus_df[:,Symbol("LE_"*string(θ)*"change")] = LE.([param_new], illus_df.age, spec = spec)
    illus_df[:,Symbol("H_"*string(θ)*"change")] = H.([param_new], illus_df.age, spec = spec)
    LE_max = maximum(illus_df[:,Symbol("LE_"*string(θ)*"change")]) + 2.5

    p1 = base_plt(illus_df, param, yscale = (0,LE_max))
    p1 = plot!(illus_df.age, illus_df[:,Symbol("μ_"*string(θ)*"change")], label = "",
        color = 1, linestyle = :dash)
    p1 = plot!(illus_df.age, illus_df[:,Symbol("S_"*string(θ)*"change")], label = "",
        color = 2, linestyle = :dash)
    p1 = plot!(illus_df.age, illus_df[:,Symbol("H_"*string(θ)*"change")], label = "",
        color = 3, linestyle = :dash)
    p1 = plot!(twinx(), illus_df[:,Symbol("LE_"*string(θ)*"change")], label = "",
        color = 4, linestyle = :dash, xlims = (-2,maximum(illus_df.age)), ylims = (0,LE_max))

    p1 = plot!(title = string(spec)*" change in "*string(θ)*" from "*
        string(round(getfield(param, θ), digits = 3))*" to "*string(round(θ_new, digits = 3)) )

    return illus_df, p1
end




"""
Function to plot the fit of model to (log) mortality data
"""
function plot_fit_year(parests, m_dist, year; log_vals = false, col = 1)

    # Extract the parameters for that year
    B = parests.mean[(parests.year .== year).*(parests.parameter .== :B)][1]
    b = parests.mean[(parests.year .== year).*(parests.parameter .== :b)][1]
    C = parests.mean[(parests.year .== year).*(parests.parameter .== :C)][1]
    c = parests.mean[(parests.year .== year).*(parests.parameter .== :c)][1]
    d = parests.mean[(parests.year .== year).*(parests.parameter .== :d)][1]
    σ = parests.mean[(parests.year .== year).*(parests.parameter .== :σ)][1]

    ages = Int.(0:length(m_dist))

    plt = plot(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
    if log_vals
        scatter!(log.(m_dist), markershape = :cross, markeralpha = 0.5,
            label = "Data ("*string(year)*")", color = :black)
        plot!(log.(siler(B,b,C,c,d, ages)),
            label = "Siler MCMC fit ("*string(year)*")", color = :blue)
        plot!(log.(siler(B,b,C,c,d, ages)) .+ 2*σ,
            label = L"\pm\ 2\ s.e.", color = :blue, linestyle = :dash)
        plot!(log.(siler(B,b,C,c,d, ages)) .- 2*σ,
            label = false, color = :blue, linestyle = :dash)
    else
        scatter!((m_dist), markershape = :cross, markeralpha = 0.5,
            label = "Data ("*string(year)*")", color = col)
        plot!((siler(B,b,C,c,d, ages)), label = "Siler MCMC fit ("*string(year)*")",
            color = col)
    end

    return plt
end


"""
Inplace version of function to plot the fit of model to (log) mortality data
"""
function plot_fit_year!(parests, m_dist, year; log_vals = false, col = 1)

    # Extract the parameters for that year
    B = parests.mean[(parests.year .== year).*(parests.parameter .== :B)][1]
    b = parests.mean[(parests.year .== year).*(parests.parameter .== :b)][1]
    C = parests.mean[(parests.year .== year).*(parests.parameter .== :C)][1]
    c = parests.mean[(parests.year .== year).*(parests.parameter .== :c)][1]
    d = parests.mean[(parests.year .== year).*(parests.parameter .== :d)][1]
    σ = parests.mean[(parests.year .== year).*(parests.parameter .== :σ)][1]

    ages = Int.(0:length(m_dist))

    plt = plot!(size = (500,300), legend = :topleft, xlab = "Age", ylab = "Mortality")
    if log_vals
        scatter!(log.(m_dist), markershape = :cross, markeralpha = 0.5,
            label = "Data ("*string(year)*")", color = :black)
        plot!(log.(siler(B,b,C,c,d, ages)),
            label = "Siler MCMC fit ("*string(year)*")", color = :blue)
        plot!(log.(siler(B,b,C,c,d, ages)) .+ 2*σ,
            label = L"\pm\ 2\ s.e.", color = :blue, linestyle = :dash)
        plot!(log.(siler(B,b,C,c,d, ages)) .- 2*σ,
            label = false, color = :blue, linestyle = :dash)
    else
        scatter!((m_dist), markershape = :cross, markeralpha = 0.5,
            label = "Data ("*string(year)*")", color = col)
        plot!((siler(B,b,C,c,d, ages)), label = "Siler MCMC fit ("*string(year)*")",
            color = col)
    end

    return plt
end


"""
Plot each Siler parameter over time
"""
function plot_siler_params(par_ests::DataFrame; forecasts = false)

    plt = plot(layout = (2,3), xrotation = 45.0, margin=3Plots.mm, size = (800,400))

    b_ests = par_ests[par_ests.parameter .== :b,:]
    plot!(b_ests.year, b_ests.median, title = L"b_{t}", label = false, color = 1, subplot = 1)
    plot!(b_ests.year, b_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(b_ests.year, b_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
    plot!(b_ests.year, b_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
    plot!(b_ests.year, b_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)
    if forecasts
        frc_ests = b_ests[b_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 1)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 1)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 1)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 1)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 1)
    end

    B_ests = par_ests[par_ests.parameter .== :B,:]
    plot!(B_ests.year, B_ests.median, title = L"B_{t}", label = false, color = 1, subplot = 2)
    plot!(B_ests.year, B_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(B_ests.year, B_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
    plot!(B_ests.year, B_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
    plot!(B_ests.year, B_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)
    if forecasts
        frc_ests = B_ests[B_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 2)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 2)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 2)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 2)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 2)
    end

    c_ests = par_ests[par_ests.parameter .== :c,:]
    plot!(c_ests.year, c_ests.median, title = L"c_{t}", label = false, color = 1, subplot = 4)
    plot!(c_ests.year, c_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(c_ests.year, c_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 4)
    plot!(c_ests.year, c_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 4)
    plot!(c_ests.year, c_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 4)
    if forecasts
        frc_ests = c_ests[c_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 4)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 4)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 4)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 4)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 4)
    end

    C_ests = par_ests[par_ests.parameter .== :C,:]
    plot!(C_ests.year, C_ests.median, title = L"C_{t}", label = false, color = 1, subplot = 5)
    plot!(C_ests.year, C_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(C_ests.year, C_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 5)
    plot!(C_ests.year, C_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 5)
    plot!(C_ests.year, C_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 5)
    if forecasts
        frc_ests = C_ests[C_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 5)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 5)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 5)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 5)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 5)
    end

    d_ests = par_ests[par_ests.parameter .== :d,:]
    plot!(d_ests.year, d_ests.median, title = L"d_{t}", label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 3)
    plot!(d_ests.year, d_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 3)
    if forecasts
        frc_ests = d_ests[d_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 3)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 3)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 3)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 3)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 3)
    end

    σ_ests = par_ests[par_ests.parameter .== :σ,:]
    plot!(σ_ests.year, σ_ests.median, title = L"\sigma_{t}", label = false, color = 1, subplot = 6)#, ylim = (0.0, 0.1))
    plot!(σ_ests.year, σ_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 6)
    plot!(σ_ests.year, σ_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 6)
    if forecasts
        frc_ests = σ_ests[σ_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 6)
        plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 6)
        plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 6)
        plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 6)
        plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 6)
    end

    return plt

end



"""
Plot the time series parameters for each Siler parameter
    model_vers defines which dynamic model (indep, justrw, i2drift, firstdiff)
"""
function plot_ts_params(par_ests::DataFrame; model_vers = :i2drift, forecasts = false)

    if model_vers == :justrw
        plt = plot(xrotation = 45.0, margin=3Plots.mm, size = (600,400))
        # Siler parameter rw variance
        σ_ests = par_ests[occursin.("σ_", string.(par_ests.parameter)),:]
        scatter!(string.(σ_ests.parameter), σ_ests.median, legend = false, color = 1)
        scatter!(string.(σ_ests.parameter), σ_ests.pc975, color = 1,markershape = :cross)
        scatter!(string.(σ_ests.parameter), σ_ests.pc025, color = 1,markershape = :cross)
        scatter!(string.(σ_ests.parameter), σ_ests.pc85, color = 1,markershape = :vline)
        scatter!(string.(σ_ests.parameter), σ_ests.pc15, color = 1,markershape = :vline)
        hline!([0], linestyle=:dot, color = :black)

    elseif model_vers == :firstdiff
        plt = plot(layout = (1,3), xrotation = 45.0, margin=3Plots.mm, size = (800,400))
        # Siler parameter rw variance
        σ_ests = par_ests[(occursin.("σ_", string.(par_ests.parameter))),:]
        scatter!(string.(σ_ests.parameter), σ_ests.median, title = "σ", legend = false, color = 1, subplot = 3)
        scatter!(string.(σ_ests.parameter), σ_ests.pc975, color = 1,markershape = :cross, subplot = 3)
        scatter!(string.(σ_ests.parameter), σ_ests.pc025, color = 1,markershape = :cross, subplot = 3)
        scatter!(string.(σ_ests.parameter), σ_ests.pc85, color = 1,markershape = :vline, subplot = 3)
        scatter!(string.(σ_ests.parameter), σ_ests.pc15, color = 1,markershape = :vline, subplot = 3)
        hline!([0], linestyle=:dot, color = :black, subplot = 3)
        # Drift  parameter
        α_ests = par_ests[(occursin.("α_", string.(par_ests.parameter))),:]
        scatter!(string.(α_ests.parameter), α_ests.median, title = "α", legend = false, color = 1, subplot = 1)
        scatter!(string.(α_ests.parameter), α_ests.pc975, color = 1,markershape = :cross, subplot = 1)
        scatter!(string.(α_ests.parameter), α_ests.pc025, color = 1,markershape = :cross, subplot = 1)
        scatter!(string.(α_ests.parameter), α_ests.pc85, color = 1,markershape = :vline, subplot = 1)
        scatter!(string.(α_ests.parameter), α_ests.pc15, color = 1,markershape = :vline, subplot = 1)
        hline!([0], linestyle=:dot, color = :black, subplot = 1)
        # Autoregressive coefficient on lagged first diff
        β_ests = par_ests[(occursin.("β_", string.(par_ests.parameter))),:]
        scatter!(string.(β_ests.parameter), β_ests.median, title = "β", legend = false, color = 1, subplot = 2)
        scatter!(string.(β_ests.parameter), β_ests.pc975, color = 1,markershape = :cross, subplot = 2)
        scatter!(string.(β_ests.parameter), β_ests.pc025, color = 1,markershape = :cross, subplot = 2)
        scatter!(string.(β_ests.parameter), β_ests.pc85, color = 1,markershape = :vline, subplot = 2)
        scatter!(string.(β_ests.parameter), β_ests.pc15, color = 1,markershape = :vline, subplot = 2)
        hline!([0], linestyle=:dot, color = :black, subplot = 2)

    elseif model_vers == :i2drift
        # Siler parameter rw variance
        σ_ests = par_ests[(occursin.("σ_", string.(par_ests.parameter)) .&
            .!occursin.("σ_α", string.(par_ests.parameter))),:]
        p_σ = scatter(string.(σ_ests.parameter), σ_ests.median, title = "σ", legend = false, color = 1)
        p_σ = scatter!(string.(σ_ests.parameter), σ_ests.pc975, color = 1,markershape = :cross)
        p_σ = scatter!(string.(σ_ests.parameter), σ_ests.pc025, color = 1,markershape = :cross)
        p_σ = scatter!(string.(σ_ests.parameter), σ_ests.pc85, color = 1,markershape = :vline)
        p_σ = scatter!(string.(σ_ests.parameter), σ_ests.pc15, color = 1,markershape = :vline)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        # Drift parameter
        ασ_ests = par_ests[(occursin.("σ_α", string.(par_ests.parameter))),:]
        p_ασ = scatter(string.(ασ_ests.parameter), ασ_ests.median, title = "σ_α", legend = false, color = 1)
        p_ασ = scatter!(string.(ασ_ests.parameter), ασ_ests.pc975, color = 1,markershape = :cross)
        p_ασ = scatter!(string.(ασ_ests.parameter), ασ_ests.pc025, color = 1,markershape = :cross)
        p_ασ = scatter!(string.(ασ_ests.parameter), ασ_ests.pc85, color = 1,markershape = :vline)
        p_ασ = scatter!(string.(ασ_ests.parameter), ασ_ests.pc15, color = 1,markershape = :vline)
        hline!([0], linestyle=:dashdot, color = :black, label = false)


        b_ests = par_ests[par_ests.parameter .== :α_b,:]
        b_plt = plot(b_ests.year, b_ests.median, label = "α_b", color = 1, xticks = false)
        plot!(b_ests.year, b_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(b_ests.year, b_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(b_ests.year, b_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(b_ests.year, b_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = b_ests[b_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        B_ests = par_ests[par_ests.parameter .== :α_B,:]
        B_plt = plot(B_ests.year, B_ests.median, label = "α_B", color = 1, xticks = false)
        plot!(B_ests.year, B_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(B_ests.year, B_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(B_ests.year, B_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(B_ests.year, B_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = B_ests[B_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        c_ests = par_ests[par_ests.parameter .== :α_c,:]
        c_plt = plot(c_ests.year, c_ests.median, label = "α_c", color = 1, xticks = false)
        plot!(c_ests.year, c_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(c_ests.year, c_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(c_ests.year, c_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(c_ests.year, c_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = c_ests[c_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        C_ests = par_ests[par_ests.parameter .== :α_C,:]
        C_plt = plot(C_ests.year, C_ests.median, label = "α_C", color = 1, xticks = false)
        plot!(C_ests.year, C_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(C_ests.year, C_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(C_ests.year, C_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(C_ests.year, C_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = C_ests[C_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        d_ests = par_ests[par_ests.parameter .== :α_d,:]
        d_plt = plot(d_ests.year, d_ests.median, label = "α_d", color = 1, xticks = false)
        plot!(d_ests.year, d_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(d_ests.year, d_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(d_ests.year, d_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(d_ests.year, d_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = d_ests[d_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        σ_ests = par_ests[par_ests.parameter .== :α_σ,:]
        σ_plt = plot(σ_ests.year, σ_ests.median, label = "α_σ", color = 1)
        plot!(σ_ests.year, σ_ests.pc975, linestyle = :dot, label = false, color = 1)
        plot!(σ_ests.year, σ_ests.pc025, linestyle = :dot, label = false, color = 1)
        plot!(σ_ests.year, σ_ests.pc85, linestyle = :dash, label = false, color = 1)
        plot!(σ_ests.year, σ_ests.pc15, linestyle = :dash, label = false, color = 1)
        hline!([0], linestyle=:dashdot, color = :black, label = false)
        if forecasts
            frc_ests = σ_ests[σ_ests.forecast .== 1,:]
            plot!(frc_ests.year, frc_ests.median, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2)
        end

        α_plts = plot(b_plt, B_plt, c_plt, C_plt, d_plt, σ_plt, layout  = grid(6,1))
        σ_plts = plot(p_σ, p_ασ, layout  = grid(2,1))
        l = @layout [a{0.5w} b{0.5w}]

        plt = plot(α_plts, σ_plts, layout = l, size = (800,800))
    end



    return plt

end


"""
Plot a decomposition of the evolution of LE or H in terms of each parameter
"""
function plot_decomp(decomp_df, variable)
  # Define variables to plot
  if variable == :LE
    b = decomp_df.Δb.*decomp_df.LE_b
    B = decomp_df.ΔB.*decomp_df.LE_B
    c = decomp_df.Δc.*decomp_df.LE_c
    C = decomp_df.ΔC.*decomp_df.LE_C
    d = decomp_df.Δd.*decomp_df.LE_d
    var = decomp_df.ΔLE_mod
  elseif variable == :H
    b = decomp_df.Δb.*decomp_df.H_b
    B = decomp_df.ΔB.*decomp_df.H_B
    c = decomp_df.Δc.*decomp_df.H_c
    C = decomp_df.ΔC.*decomp_df.H_C
    d = decomp_df.Δd.*decomp_df.H_d
    var = decomp_df.ΔH_mod
  end
  # Plot decomposition
  p1 = groupedbar(decomp_df.year, [b B c C d], label=["b" "B" "c" "C" "d"],
    bar_position = :stack, linecolor=nothing,
    xlabel = "Year", ylabel = "Δ"*string(variable))
  p1 = plot!(decomp_df.year, var, color = :black, label = false)
  p1 = hline!([0,0], color = :black, linestyle = :dash, label = false)

  return p1
end


"""
Function to plot model implied LE and H
"""
function plot_LE_H(par_ests; forecasts = false, bands = false)
    plt = plot(layout = (1,2), xrotation = 45.0, margin=3Plots.mm, size = (800,400))

    LE_ests = par_ests[par_ests.parameter .== :LE,:]
    plot!(LE_ests.year, LE_ests.median, title = L"LE_{t}", label = false, color = 1, subplot = 1)
    if bands
        plot!(LE_ests.year, LE_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 1)
        plot!(LE_ests.year, LE_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 1)
        plot!(LE_ests.year, LE_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 1)
        plot!(LE_ests.year, LE_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 1)
    end
    if forecasts
        frc_ests = LE_ests[LE_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 1)
        if bands
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 1)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 1)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 1)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 1)
        end
    end

    H_ests = par_ests[par_ests.parameter .== :H,:]
    plot!(H_ests.year, H_ests.median, title = L"H_{t}", label = false, color = 1, subplot = 2)
    if bands
        plot!(H_ests.year, H_ests.pc975, linestyle = :dot, label = false, color = 1, subplot = 2)
        plot!(H_ests.year, H_ests.pc025, linestyle = :dot, label = false, color = 1, subplot = 2)
        plot!(H_ests.year, H_ests.pc85, linestyle = :dash, label = false, color = 1, subplot = 2)
        plot!(H_ests.year, H_ests.pc15, linestyle = :dash, label = false, color = 1, subplot = 2)
    end
    if forecasts
        frc_ests = H_ests[H_ests.forecast .== 1,:]
        plot!(frc_ests.year, frc_ests.median, label = false, color = 2, subplot = 2)
        if bands
            plot!(frc_ests.year, frc_ests.pc975, linestyle = :dot, label = false, color = 2, subplot = 2)
            plot!(frc_ests.year, frc_ests.pc025, linestyle = :dot, label = false, color = 2, subplot = 2)
            plot!(frc_ests.year, frc_ests.pc85, linestyle = :dash, label = false, color = 2, subplot = 2)
            plot!(frc_ests.year, frc_ests.pc15, linestyle = :dash, label = false, color = 2, subplot = 2)
        end
    end

    return plt

end
