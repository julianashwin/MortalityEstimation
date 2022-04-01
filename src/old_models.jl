"""
Dynamic Siler model with random walk in log parameters
"""
@model function log_siler_justrw(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Log parameters for random walk model
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.1),6)
    # First period from priors
    lB[1] ~ Normal(log(5), 2.0)
    lb[1] ~ Normal(log(1), 1.0)
    lC[1] ~ Normal(log(10), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.025), 1.0)
    lσ[1] ~ Normal(log(0.1), 1.0)
    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*ages .- exp(lB[1])) .+
        exp.(exp(lc[1]).*ages .- exp(lC[1])) .+ exp(ld[1])
    lm_vars = exp(lσ[1]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ)
    # Loop through random walk process
    for tt in 2:T
        # Calculate updated variances
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
        μs = exp.(-exp(lb[tt]).*ages .- exp(lB[tt])) .+
            exp.(exp(lc[tt]).*ages .- exp(lC[tt])) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ)
    end
end




"""
Dynamic Siler model with parameters evolving as first differences
    - drift term
    - Non-diagonal covariance matrix (not working)
    - AR(1) with trend (not working)

Note that mean of InverseWishart is S/(m - p - 1) where S is scale matrix, m is d.o.f., p is dimensions
"""
@model function log_siler_dyn_firstdiff(lm_data, ages)
    # Dimensions
    T::Int64 = length(lm_data)
    N::Int64 = length(ages)
    # Parameters
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Priors on variance terms for parameter time series
    σ_pars ~ filldist(InverseGamma(2, 0.05),6)
    # Priors on constant
    α_pars ~ filldist(Normal(0., 0.05), 6)
    # Autoregressive coefficient
    β_pars ~ filldist(Normal(0.0, 0.3), 6)

    # Independent draws for the first two periods
    for tt in 1:2
        lB[tt] ~ Normal(log(10), 2.0)
        lb[tt] ~ Normal(log(2), 1.0)
        lC[tt] ~ Normal(log(10), 2.0)
        lc[tt] ~ Normal(log(0.1), 1.0)
        ld[tt] ~ Normal(log(0.01), 1.0)
        lσ[tt] ~ Normal(log(0.01), 1.0)
        # Find mean using the siler mortality function
        μs = exp.(-(exp(lb[tt]).*ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*ages .- exp(lC[tt])) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)
    end

    # Loop through random walk process
    for tt in 3:T
        # Calculate updated parameter means
        μ_B = lB[tt-1] + α_pars[1] + β_pars[1]*(lB[tt-1] - lB[tt-2])
        μ_b = lb[tt-1] + α_pars[2] + β_pars[2]*(lb[tt-1] - lb[tt-2])
        μ_C = lC[tt-1] + α_pars[3] + β_pars[3]*(lC[tt-1] - lC[tt-2])
        μ_c = lc[tt-1] + α_pars[4] + β_pars[4]*(lc[tt-1] - lc[tt-2])
        μ_d = ld[tt-1] + α_pars[5] + β_pars[5]*(ld[tt-1] - ld[tt-2])
        μ_σ = lσ[tt-1] + α_pars[6] + β_pars[6]*(lσ[tt-1] - lσ[tt-2])
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
        μs = exp.(-(exp(lb[tt]).*ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*ages .- exp(lC[tt])) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)
    end
end
