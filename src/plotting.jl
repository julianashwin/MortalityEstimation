
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
