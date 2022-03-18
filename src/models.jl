"""
Turing models for Siler function mortality
"""



"""
Static Siler model on non-logged data
"""
@model function siler_static(m_dist, ages)
    # The number of observations.
    N = length(m_dist)
    # Our prior beliefs
    B ~ LogNormal(log(10), 2.0)
    b ~ LogNormal(log(2), 1.0)
    C ~ LogNormal(log(10), 2.0)
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
    μs = exp.(- exp(lb).* ages .+ exp(lB)) .+ exp.(exp(lc) .* ages .- exp(lC)) .+ exp(ld)
    m_var = exp(σ)
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
    μs = exp.(- exp(lb).* ages .+ exp(lB)) .+ exp.(exp(lc) .* ages .- exp(lC)) .+ exp(ld)
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
    C ~ filldist(LogNormal(log(10), 2.0), T)
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
        μs = exp.(-exp(lb[tt]).*ages .+ exp(lB[tt])) .+
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
    lC[1] ~ Normal(log(10), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.025), 1.0)
    lσ[1] ~ Normal(log(0.1), 1.0)
    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*ages .+ exp(lB[1])) .+
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
        μs = exp.(-exp(lb[tt]).*ages .+ exp(lB[tt])) .+
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
    σ_par ~ filldist(InverseGamma(2, 0.05),6)
    # Correlation matrix of shocks to parameter time series
    #ρ_block ~ LKJ(6, 4)
    # Prior on variance terms for parameter time series
    #Σ_ϵ = Diagonal(σ_par)#quad_form_diag(ρ_block, σ_par)
    # Priors on constant
    α_pars ~ filldist(Normal(0., 0.05), 6)
    # Priors on time trend terms
    #τ_pars ~ filldist(Normal(0, 0.05), 6)
    # Autoregressive coefficient
    β_pars ~ filldist(Normal(0.0, 0.3), 6)

    # Independent draws for the first two periods
    for tt in 1:2
        lB[tt] ~ Normal(log(10), 2.0)
        lb[tt] ~ Normal(log(2), 1.0)
        lC[tt] ~ Normal(log(10), 2.0)
        lc[tt] ~ Normal(log(0.1), 1.0)
        ld[tt] ~ Normal(log(0.025), 1.0)
        lσ[tt] ~ Normal(log(0.1), 1.0)
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*ages .+ exp(lB[tt])) .+
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
        var_B = max(σ_par[1], 1e-8)
        var_b = max(σ_par[2], 1e-8)
        var_C = max(σ_par[3], 1e-8)
        var_c = max(σ_par[4], 1e-8)
        var_d = max(σ_par[5], 1e-8)
        var_σ = max(σ_par[6], 1e-8)
        # Update parameters
        lB[tt] ~ Normal(μ_B, var_B)
        lb[tt] ~ Normal(μ_b, var_b)
        lC[tt] ~ Normal(μ_C, var_C)
        lc[tt] ~ Normal(μ_c, var_c)
        ld[tt] ~ Normal(μ_d, var_d)
        lσ[tt] ~ Normal(μ_σ, var_σ)

        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*ages .+ exp(lB[tt])) .+
            exp.(exp(lc[tt]).*ages .- exp(lC[tt])) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)
    end
end









"""
Dynamic Siler model with extensions
    - time-varying drift term
    - Non-diagonal covariance matrix (not working)
    - AR(1) with trend (not working)

Note that mean of InverseWishart is S/(m - p - 1) where S is scale matrix, m is d.o.f., p is dimensions
"""
@model function log_siler_dyn_ext(lm_data, ages)
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
    σ_par ~ filldist(InverseGamma(2, 0.05),6)
    # Correlation matrix of shocks to parameter time series
    # Priors on constant
    α_B = Vector(undef, T-1)
    α_b = Vector(undef, T-1)
    α_C = Vector(undef, T-1)
    α_c = Vector(undef, T-1)
    α_d = Vector(undef, T-1)
    α_σ = Vector(undef, T-1)
    # Variance of drift random walk
    σ_αB ~ InverseGamma(2, 0.01)
    σ_αb ~ InverseGamma(2, 0.01)
    σ_αC ~ InverseGamma(2, 0.01)
    σ_αc ~ InverseGamma(2, 0.01)
    σ_αd ~ InverseGamma(2, 0.01)
    σ_ασ ~ InverseGamma(2, 0.01)
    # Priors on time trend terms
    #τ_pars ~ filldist(Normal(0, 0.05), 6)
    # Autoregressive coefficient
    #β_pars ~ filldist(Uniform(0.0, 0.99), 6
    # First period from priors
    lB[1] ~ Normal(log(10), 2.0)
    lb[1] ~ Normal(log(2), 1.0)
    lC[1] ~ Normal(log(120), 2.0)
    lc[1] ~ Normal(log(0.1), 1.0)
    ld[1] ~ Normal(log(0.025), 1.0)
    lσ[1] ~ Normal(log(0.1), 1.0)
    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*ages .+ exp(lB[1])) .+
        exp.(exp(lc[1]).*ages .- exp(lC[1])) .+ exp(ld[1])
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
        #μ_B = α_pars[1] + β_pars[1]*lB[tt-1] + τ_pars[1]*(tt-1)
        #μ_b = α_pars[2] + β_pars[2]*lb[tt-1] + τ_pars[2]*(tt-1)
        #μ_C = α_pars[3] + β_pars[3]*lC[tt-1] + τ_pars[3]*(tt-1)
        #μ_c = α_pars[4] + β_pars[4]*lc[tt-1] + τ_pars[4]*(tt-1)
        #μ_d = α_pars[5] + β_pars[5]*ld[tt-1] + τ_pars[5]*(tt-1)
        #μ_σ = α_pars[6] + β_pars[6]*lσ[tt-1] + τ_pars[6]*(tt-1)
        μ_B = α_B[tt-1] + lB[tt-1]
        μ_b = α_b[tt-1] + lb[tt-1]
        μ_C = α_C[tt-1] + lC[tt-1]
        μ_c = α_c[tt-1] + lc[tt-1]
        μ_d = α_d[tt-1] + ld[tt-1]
        μ_σ = α_σ[tt-1] + lσ[tt-1]
        # Updated variance parameters for the mortality curve
        var_B = max(σ_par[1], 1e-8)
        var_b = max(σ_par[2], 1e-8)
        var_C = max(σ_par[3], 1e-8)
        var_c = max(σ_par[4], 1e-8)
        var_d = max(σ_par[5], 1e-8)
        var_σ = max(σ_par[6], 1e-8)
        # Update parameters
        lB[tt] ~ Normal(μ_B, var_B)
        lb[tt] ~ Normal(μ_b, var_b)
        lC[tt] ~ Normal(μ_C, var_C)
        lc[tt] ~ Normal(μ_c, var_c)
        ld[tt] ~ Normal(μ_d, var_d)
        lσ[tt] ~ Normal(μ_σ, var_σ)

        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*ages .+ exp(lB[tt])) .+
            exp.(exp(lc[tt]).*ages .- exp(lC[tt])) .+ exp(ld[tt])
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
# Priors on variance terms for parameter time series
σ_B ~ InverseGamma(2, 0.3)
σ_b ~ InverseGamma(2, 0.3)
σ_C ~ InverseGamma(2, 0.3)
σ_c ~ InverseGamma(2, 0.3)
σ_d ~ InverseGamma(2, 0.3)
σ_σ ~ InverseGamma(2, 0.3)

# Priors on covariance terms for parameter time series
τ_Bb ~ Uniform(-0.99, 0.99)
τ_BC ~ Uniform(-0.99, 0.99)
τ_Bc ~ Uniform(-0.99, 0.99)
τ_Bd ~ Uniform(-0.99, 0.99)
τ_Bσ ~ Uniform(-0.99, 0.99)
τ_bC ~ Uniform(-0.99, 0.99)
τ_bc ~ Uniform(-0.99, 0.99)
τ_bd ~ Uniform(-0.99, 0.99)
τ_bσ ~ Uniform(-0.99, 0.99)
τ_Cc ~ Uniform(-0.99, 0.99)
τ_Cd ~ Uniform(-0.99, 0.99)
τ_Cσ ~ Uniform(-0.99, 0.99)
τ_cd ~ Uniform(-0.99, 0.99)
τ_cσ ~ Uniform(-0.99, 0.99)
τ_dσ ~ Uniform(-0.99, 0.99)

Σ_ϵ = Symmetric(
      [σ_B^2 + 1e-6    τ_Bb*σ_B*σ_b    τ_BC*σ_B*σ_C    τ_Bc*σ_B*σ_c    τ_Bd*σ_B*σ_d    τ_Bσ*σ_B*σ_σ ;
       τ_Bb*σ_B*σ_b    σ_b^2 + 1e-6    τ_bC*σ_b*σ_C    τ_bc*σ_b*σ_c    τ_bd*σ_b*σ_d    τ_bσ*σ_b*σ_σ ;
       τ_BC*σ_B*σ_C    τ_bC*σ_b*σ_C    σ_C^2 + 1e-6    τ_Cc*σ_C*σ_c    τ_Cd*σ_C*σ_d    τ_Cσ*σ_C*σ_σ ;
       τ_Bc*σ_B*σ_c    τ_bc*σ_b*σ_c    τ_Cc*σ_C*σ_c    σ_c^2 + 1e-6    τ_cd*σ_c*σ_d    τ_cσ*σ_c*σ_σ ;
       τ_Bd*σ_B*σ_d    τ_bd*σ_b*σ_d    τ_Cd*σ_C*σ_d    τ_cd*σ_c*σ_d    σ_d^2 + 1e-6    τ_dσ*σ_d*σ_σ;
       τ_Bσ*σ_B*σ_σ    τ_bσ*σ_b*σ_σ    τ_Cσ*σ_C*σ_σ    τ_cσ*σ_c*σ_σ    τ_dσ*σ_d*σ_σ    σ_σ^2 + 1e-6])

"""
