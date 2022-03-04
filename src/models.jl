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





"""
Dynamic Siler model with extensions
    - Non-diagonal covariance matrix
"""
@model function log_siler_dyn_ext(lm_data, ages)
    # Dimensions
    T = length(lm_data)
    N = length(ages)
    # Parameters
    lB = Vector(undef, T)
    lb = Vector(undef, T)
    lC = Vector(undef, T)
    lc = Vector(undef, T)
    ld = Vector(undef, T)
    lσ = Vector(undef, T)
    # Matrix of all parameters (will be drawn from MVNormal)
    pars = Matrix(undef,6,T)
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
           [σ_B^2    τ_Bb*σ_B*σ_b    τ_BC*σ_B*σ_C    τ_Bc*σ_B*σ_c    τ_Bd*σ_B*σ_d    τ_Bσ*σ_B*σ_σ ;
           τ_Bb*σ_B*σ_b    σ_b^2    τ_bC*σ_b*σ_C    τ_bc*σ_b*σ_c    τ_bd*σ_b*σ_d    τ_bσ*σ_b*σ_σ ;
           τ_BC*σ_B*σ_C    τ_bC*σ_b*σ_C    σ_C^2    τ_Cc*σ_C*σ_c    τ_Cd*σ_C*σ_d    τ_Cσ*σ_C*σ_σ ;
           τ_Bc*σ_B*σ_c    τ_bc*σ_b*σ_c    τ_Cc*σ_C*σ_c    σ_c^2    τ_cd*σ_c*σ_d    τ_cσ*σ_c*σ_σ ;
           τ_Bd*σ_B*σ_d    τ_bd*σ_b*σ_d    τ_Cd*σ_C*σ_d    τ_cd*σ_c*σ_d    σ_d^2    τ_dσ*σ_d*σ_σ;
           τ_Bσ*σ_B*σ_σ    τ_bσ*σ_b*σ_σ    τ_Cσ*σ_C*σ_σ    τ_cσ*σ_c*σ_σ    τ_dσ*σ_d*σ_σ    σ_σ^2 ])



    # Priors on drift terms
    #μ_pars ~ filldist(Normal(0, 0.1), 6)
    # First period from priors
    μ_par = log.([10.;2.;120.;0.1;0.25;0.1])
    pars[:,1] ~ MvNormal(μ_par, Diagonal([2.;1.;2.;1.;1.;1.]))
    lB[1] = pars[1,1]
    lb[1] = pars[2,1]
    lC[1] = pars[3,1]
    lc[1] = pars[4,1]
    ld[1] = pars[5,1]
    lσ[1] = pars[6,1]

    # Find mean using the siler mortality function
    μs = exp.(-exp(lb[1]).*(ages .+ exp(lB[1]))) .+
        exp.(exp(lc[1]).*(ages .- exp(lC[1]))) .+ exp(ld[1])
    lm_vars = exp(lσ[1]).*ones(N)
    lm_vars[lm_vars.<= 1e-10] .= 1e-10
    # Variance matrix
    Σ_m = Diagonal(lm_vars)
    # Draw from truncated normal dist
    lm_data[1] ~ MvNormal(log.(μs), Σ_m)
    # Loop through random walk process
    for tt in 2:T
        # Calculate updated means
        μ_par = [lB[tt-1];
                 lb[tt-1];
                 lC[tt-1];
                 lc[tt-1];
                 ld[tt-1];
                 lσ[tt-1]]
        # Update parameters
        pars[:,tt] ~ MvNormal(μ_par, Σ_ϵ)
        lB[tt] = pars[1,tt]
        lb[tt] = pars[2,tt]
        lC[tt] = pars[3,tt]
        lc[tt] = pars[4,tt]
        ld[tt] = pars[5,tt]
        lσ[tt] = pars[6,tt]
        # Find mean using the siler mortality function
        μs = exp.(-exp(lb[tt]).*(ages .+ exp(lB[tt]))) .+
            exp.(exp(lc[tt]).*(ages .- exp(lC[tt]))) .+ exp(ld[tt])
        lm_vars = exp(lσ[tt]).*ones(N)
        lm_vars[lm_vars.<= 1e-10] .= 1e-10
        # Variance matrix
        Σ_m = Diagonal(lm_vars)
        # Draw from truncated normal dist
        lm_data[tt] ~ MvNormal(log.(μs), Σ_m)
    end
end
