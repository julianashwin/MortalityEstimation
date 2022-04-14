"""
Turing models for Siler function mortality
    In all cases, set up the model as in Colchero: μ = e^{-ba -B} + e^{ca - C} + d
    We can then transform the results using options in the extract_variables function (helpers.jl)
"""



"""
Static Siler model on non-logged data
"""
@model function siler_static(m_dist, ages)
    # The number of observations.
    N = length(m_dist)
    # Our prior beliefs
    B ~ LogNormal(log(2), 2.0)
    b ~ LogNormal(log(1), 1.0)
    C ~ LogNormal(log(80), 2.0)
    c ~ LogNormal(log(0.1), 1.0)
    d ~ LogNormal(log(0.01), 1.0)
    σ ~ LogNormal(log(0.001), 1.0)
    # Find mean using the siler mortality function
    μs = b.*exp.(- b.*(ages .+ B)) .+ c.*exp.(c.*(ages .- C)) .+ d
    m_var = σ
    #m_vars[m_vars.<= 1e-10] .= 1e-10
    # Draw from normal dist
    for nn in 1:N
        m_dist[nn] ~ LogNormal(log(μs[nn]), m_var)
    end
    #m_dist = exp.(lm_dist)
end



"""
Static Siler model on log data
"""
@model function log_siler_static(lm_dist, ages)
    # The number of observations.
    N = length(lm_dist)
    # Our prior beliefs
    B ~ LogNormal(log(2), 2.0)
    b ~ LogNormal(log(1), 1.0)
    C ~ LogNormal(log(80), 2.0)
    c ~ LogNormal(log(0.1), 1.0)
    d ~ LogNormal(log(0.01), 1.0)
    σ ~ LogNormal(log(0.001), 1.0)
    # Find mean using the siler mortality function
    μs = b.*exp.(- b.* (ages .+ B)) .+ c.*exp.(c .* (ages .- C)) .+ d
    m_vars = σ.*ones(N)
    #m_vars[m_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ = Diagonal(m_vars)
    # Draw from normal dist
    lm_dist ~ MvNormal(log.(μs), Σ)
    #m_dist = exp.(lm_dist)
end





"""
Multiple independent Siler models
"""
@model function log_siler_indep(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Our prior beliefs
    B ~ filldist(LogNormal(log(2), 2.0), T)
    b ~ filldist(LogNormal(log(1), 1.0), T)
    C ~ filldist(LogNormal(log(80), 2.0), T)
    c ~ filldist(LogNormal(log(0.1), 1.0), T)
    d ~ filldist(LogNormal(log(0.01), 1.0), T)
    σ ~ filldist(LogNormal(log(0.001), 1.0), T)
    for tt in 1:T
        # Find mean using the siler mortality function
        μs = b[tt].*exp.(-b[tt].*(ages .+ B[tt])) .+
            c[tt].*exp.(c[tt].*(ages .- C[tt])) .+ d[tt]
        lm_vars = σ[tt].*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end




"""
Dynamic Siler model with time-varying drift
    - time-varying drift term
    - Non-diagonal covariance matrix (not working)
"""
@model function log_siler_dyn_i2drift(lm_data, ages)
    # Dimensions
    T::Int64 = length(lm_data)
    N::Int64 = length(ages)
    # Siler parameters
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.05),6)
    # Time varying drift
    α_B = Vector(undef, T-1)
    α_b = Vector(undef, T-1)
    α_C = Vector(undef, T-1)
    α_c = Vector(undef, T-1)
    α_d = Vector(undef, T-1)
    α_σ = Vector(undef, T-1)
    # Variance of drift random walk
    σ_αB ~ InverseGamma(2, 0.02)
    σ_αb ~ InverseGamma(2, 0.02)
    σ_αC ~ InverseGamma(2, 0.02)
    σ_αc ~ InverseGamma(2, 0.02)
    σ_αd ~ InverseGamma(2, 0.02)
    σ_ασ ~ InverseGamma(2, 0.02)
    # First period from priors
    lB[1] ~ Normal(log(5), 2.0)
    lb[1] ~ Normal(log(1), 1.0)
    lC[1] ~ Normal(log(80), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.01), 1.0)
    lσ[1] ~ Normal(log(0.01), 1.0)
    # Find mean using the siler mortality function
    μs = exp(lb[1]).*exp.(-exp(lb[1]).*(ages .+ exp(lB[1]))) .+
        exp(lc[1]).*exp.(exp(lc[1]).*(ages .- exp(lC[1]))) .+ exp(ld[1])
    lm_vars = exp(lσ[1]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ_m = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ_m)

    # First period drift terms
    #β_pars ~ filldist(Uniform(0.0, 0.99), 6
    α_B[1] ~ Normal(0., 0.02)
    α_b[1] ~ Normal(0., 0.02)
    α_C[1] ~ Normal(0., 0.02)
    α_c[1] ~ Normal(0., 0.02)
    α_d[1] ~ Normal(0., 0.02)
    α_σ[1] ~ Normal(0., 0.02)

    # Loop through random walk process
    for tt in 2:T
        # Calculate updated parameter means
        μ_B = α_B[tt-1] + lB[tt-1]
        μ_b = α_b[tt-1] + lb[tt-1]
        μ_C = α_C[tt-1] + lC[tt-1]
        μ_c = α_c[tt-1] + lc[tt-1]
        μ_d = α_d[tt-1] + ld[tt-1]
        μ_σ = α_σ[tt-1] + lσ[tt-1]
        # Updated variance parameters for the mortality curve
        var_B = max(σ_pars[1], 1e-8)
        var_b = max(σ_pars[2], 1e-8)
        var_C = max(σ_pars[3], 1e-8)
        var_c = max(σ_pars[4], 1e-8)
        var_d = max(σ_pars[5], 1e-8)
        var_σ = max(σ_pars[6], 1e-8)
        # Update parameters
        lB[tt] ~ Normal(μ_B, var_B)
        lb[tt] ~ Normal(μ_b, var_b)
        lC[tt] ~ Normal(μ_C, var_C)
        lc[tt] ~ Normal(μ_c, var_c)
        ld[tt] ~ Normal(μ_d, var_d)
        lσ[tt] ~ Normal(μ_σ, var_σ)
        # Find mean using the siler mortality function
        μs = exp(lb[tt]).*exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp(lc[tt]).*exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)

        # Update the r.w. with drift term
        if tt < T
            α_B[tt] ~ Normal(α_B[tt-1], σ_αB)
            α_b[tt] ~ Normal(α_b[tt-1], σ_αb)
            α_C[tt] ~ Normal(α_C[tt-1], σ_αC)
            α_c[tt] ~ Normal(α_c[tt-1], σ_αc)
            α_d[tt] ~ Normal(α_d[tt-1], σ_αd)
            α_σ[tt] ~ Normal(α_σ[tt-1], σ_ασ)
        end
    end
end









"""
Dynamic Siler model with time-varying drift
    - time-varying drift term
    - Non-diagonal covariance matrix (not working)
"""
@model function log_siler_dyn_i2drift_cov(lm_data, ages)
    # Dimensions
    T::Int64 = length(lm_data)
    N::Int64 = length(ages)
    # Siler parameters
    lpars = [Vector(undef,6) for _ in 1:T]
    #lB = Vector(undef, T)
    #lb = Vector(undef, T)
    #lC = Vector(undef, T)
    #lc = Vector(undef, T)
    #ld = Vector(undef, T)
    #lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.05),6)
    ρ_bB ~ Uniform(-0.5,0.5)
    ρ_cC ~ Uniform(-0.5,0.5)
    ρ_Cd ~ Uniform(-0.5,0.5)
    ρ_cd ~ Uniform(-0.5,0.5)
    #ρ_Cd = 0.6
    #ρ_cd = -0.6
    #ρ_Cc = 0.6
    Σ_ϵ = Matrix(Diagonal(σ_pars))
    Σ_ϵ[1,2] = ρ_bB*sqrt(σ_pars[1])*sqrt(σ_pars[2])
    Σ_ϵ[2,1] = ρ_bB*sqrt(σ_pars[1])*sqrt(σ_pars[2])
    Σ_ϵ[3,4] = ρ_cC*sqrt(σ_pars[3])*sqrt(σ_pars[4])
    Σ_ϵ[4,3] = ρ_cC*sqrt(σ_pars[3])*sqrt(σ_pars[4])
    Σ_ϵ[3,5] = ρ_Cd*sqrt(σ_pars[3])*sqrt(σ_pars[5])
    Σ_ϵ[5,3] = ρ_Cd*sqrt(σ_pars[3])*sqrt(σ_pars[5])
    Σ_ϵ[4,5] = ρ_cd*sqrt(σ_pars[4])*sqrt(σ_pars[5])
    Σ_ϵ[5,4] = ρ_cd*sqrt(σ_pars[4])*sqrt(σ_pars[5])
    #eigen(Σ_ϵ)
    #PDMat(Σ_ϵ)
    # Time varying drift
    α_B = Vector(undef, T-1)
    α_b = Vector(undef, T-1)
    α_C = Vector(undef, T-1)
    α_c = Vector(undef, T-1)
    α_d = Vector(undef, T-1)
    α_σ = Vector(undef, T-1)
    # Variance of drift random walk
    σ_αB ~ InverseGamma(2, 0.02)
    σ_αb ~ InverseGamma(2, 0.02)
    σ_αC ~ InverseGamma(2, 0.02)
    σ_αc ~ InverseGamma(2, 0.02)
    σ_αd ~ InverseGamma(2, 0.02)
    σ_ασ ~ InverseGamma(2, 0.02)
    # First period from priors
    lpars[1][1] ~ Normal(log(5), 2.0)
    lpars[1][2] ~ Normal(log(1), 1.0)
    lpars[1][3] ~ Normal(log(80), 2.0)
    lpars[1][4] ~ Normal(log(0.1), 1.0)
    lpars[1][5] ~ Normal(log(0.01), 1.0)
    lpars[1][6] ~ Normal(log(0.01), 1.0)
    #lB[1] ~ Normal(log(5), 2.0)
    #lb[1] ~ Normal(log(1), 1.0)
    #lC[1] ~ Normal(log(80), 2.0)
    #lc[1] ~ Normal(log(0.1), 1.0)
    #ld[1] ~ Normal(log(0.01), 1.0)
    #lσ[1] ~ Normal(log(0.05), 1.0)
    # Find mean using the siler mortality function
    μs = exp(lpars[1][2]).*exp.(-exp(lpars[1][1]).*(ages .+ exp(lpars[1][2]))) .+
        exp(lpars[1][4]).*exp.(exp(lpars[1][4]).*(ages .- exp(lpars[1][3]))) .+ exp(lpars[1][5])
    lm_vars = exp(lpars[1][6]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ_m = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ_m)

    # First period drift terms
    #β_pars ~ filldist(Uniform(0.0, 0.99), 6
    α_B[1] ~ Normal(0., 0.02)
    α_b[1] ~ Normal(0., 0.02)
    α_C[1] ~ Normal(0., 0.02)
    α_c[1] ~ Normal(0., 0.02)
    α_d[1] ~ Normal(0., 0.02)
    α_σ[1] ~ Normal(0., 0.02)

    # Loop through random walk process
    for tt in 2:T
        # Calculate updated parameter means
        αs = [α_B[tt-1], α_b[tt-1], α_C[tt-1], α_c[tt-1], α_d[tt-1], α_σ[tt-1]]
        μ_pars = αs .+ lpars[tt-1]
        # Updated variance parameters for the mortality curve
        #var_B = max(σ_pars[1], 1e-8)
        #var_b = max(σ_pars[2], 1e-8)
        #var_C = max(σ_pars[3], 1e-8)
        #var_c = max(σ_pars[4], 1e-8)
        #var_d = max(σ_pars[5], 1e-8)
        #var_σ = max(σ_pars[6], 1e-8)
        # Update parameters
        lpars[tt] ~ MvNormal(μ_pars, Σ_ϵ)
        #lB[tt] ~ Normal(μ_B, var_B)
        #lb[tt] ~ Normal(μ_b, var_b)
        #lC[tt] ~ Normal(μ_C, var_C)
        #lc[tt] ~ Normal(μ_c, var_c)
        #ld[tt] ~ Normal(μ_d, var_d)
        #lσ[tt] ~ Normal(μ_σ, var_σ)
        # Find mean using the siler mortality function
        μs = exp(lpars[tt][2]).*exp.(-exp(lpars[tt][1]).*(ages .+ exp(lpars[tt][2]))) .+
            exp(lpars[tt][4]).*exp.(exp(lpars[tt][4]).*(ages .- exp(lpars[tt][3]))) .+ exp(lpars[tt][5])
        lm_vars = exp(lpars[tt][6]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)

        # Update the r.w. with drift term
        if tt < T
            α_B[tt] ~ Normal(α_B[tt-1], σ_αB)
            α_b[tt] ~ Normal(α_b[tt-1], σ_αb)
            α_C[tt] ~ Normal(α_C[tt-1], σ_αC)
            α_c[tt] ~ Normal(α_c[tt-1], σ_αc)
            α_d[tt] ~ Normal(α_d[tt-1], σ_αd)
            α_σ[tt] ~ Normal(α_σ[tt-1], σ_ασ)
        end
    end
end
