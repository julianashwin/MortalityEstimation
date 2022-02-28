"""
Turing models for Siler function mortality
"""


"""
Static Siler model
"""
@model function log_siler_static(lm_dist, ages)
    # The number of observations.
    N = length(lm_dist)
    # Our prior beliefs
    B ~ LogNormal(log(10), 2.0)
    b ~ LogNormal(log(2), 1.0)
    C ~ LogNormal(log(120), 2.0)
    c ~ LogNormal(log(0.1), 1.0)
    d ~ LogNormal(log(0.025), 1.0)
    σ ~ LogNormal(log(0.001), 1.0)
    # Define the logged parameters, which should be normally distributed
    lB = log(B)
    lb = log(b)
    lC = log(C)
    lc = log(c)
    ld = log(d)
    lσ = log(σ)
    # Find mean using the siler mortality function
    μs = exp.(- exp(lb).* (ages .+ exp(lB))) .+ exp.(exp(lc) .* (ages .- exp(lC))) .+ exp(ld)
    m_vars = exp(σ).*ones(N)
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
## Define Turing model
@model function log_siler_indep(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Our prior beliefs
    B ~ filldist(LogNormal(log(10), 2.0), T)
    b ~ filldist(LogNormal(log(2), 1.0), T)
    C ~ filldist(LogNormal(log(120), 2.0), T)
    c ~ filldist(LogNormal(log(0.1), 1.0), T)
    d ~ filldist(LogNormal(log(0.025), 1.0), T)
    σ ~ filldist(LogNormal(log(0.1), 1.0), T)
    # Define the logged parameters, which should be normally distributed
    lB = log.(B)
    lb = log.(b)
    lC = log.(C)
    lc = log.(c)
    ld = log.(d)
    lσ = log.(σ)

    for tt in 1:T
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end






"""
Dynamic Siler model with random walk in log parameters
"""
@model function log_siler_dyn(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Parameters
    B = Vector(undef, T)
    b = Vector(undef, T)
    C = Vector(undef, T)
    c = Vector(undef, T)
    d = Vector(undef, T)
    σ = Vector(undef, T)
    # Logged parameters for random walk model
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.1),6)
    # Priors on drift terms
    #μ_pars ~ filldist(Normal(0, 0.1), 6)

    # First period from priors
    lB[1] ~ Normal(log(10), 2.0)
    lb[1] ~ Normal(log(2), 1.0)
    lC[1] ~ Normal(log(120), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.025), 1.0)
    lσ[1] ~ Normal(log(0.1), 1.0)
    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*(ages .+ exp(lB[1]))) .+
        exp.(exp(lc[1]).*(ages .- exp(lC[1]))) .+ exp(ld[1])
    lm_vars = exp(lσ[1]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ)
    # Loop through random walk process
    for tt in 2:T
        # Calculate updated variances
        #Σ_pars = Diagonal(σ_pars)
        #lmean_pars = [lB[tt-1], lb[tt-1], lC[tt-1], lc[tt-1], ld[tt-1], lσ[tt-1]]
        var_B = max(σ_pars[1], 1e-8)
        var_b = max(σ_pars[2], 1e-8)
        var_C = max(σ_pars[3], 1e-8)
        var_c = max(σ_pars[4], 1e-8)
        var_d = max(σ_pars[5], 1e-8)
        var_σ = max(σ_pars[6], 1e-8)
        # Update parameters
        lB[tt] ~ Normal(lB[tt-1], var_B)
        lb[tt] ~ Normal(lb[tt-1], var_b)
        lC[tt] ~ Normal(lC[tt-1], var_C)
        lc[tt] ~ Normal(lc[tt-1], var_c)
        ld[tt] ~ Normal(ld[tt-1], var_d)
        lσ[tt] ~ Normal(lσ[tt-1], var_σ)
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end
