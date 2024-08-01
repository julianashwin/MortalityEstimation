
"""
Basic Siler mortality function
"""
function siler(param::SilerParam, age::Real; spec::Symbol = :Bergeron )
    @unpack b, B, c, C, d = param
	if spec == :Colchero
    	μ = exp.(- (b.* age .+ B)) .+ exp.(c .* age.- C) .+ d
	elseif spec == :Scott
		μ = exp.(- b.* (age .+ B)) .+ exp.(c .* (age.- C)) .+ d
	elseif spec == :Bergeron
		μ = b.*exp.(- b.* (age + B)) .+ c.*exp.(c .* (age.- C)) .+ d
	end
	μ = min.(μ, 1.0)

    return μ
end

function siler_thresh(param::SilerParam, age::Real, ā::Int, g::Real)
    @unpack b, B, c, C, d = param
	if age < ā
		μ = b.*exp.(- b.* (age + B)) .+ c.*exp.(c .* (age.- C)) .+ d
	else
		μ = g.*age
	end

	μ = min.(μ, 1.0)

    return μ
end




"""
Basic Heligman Pollard mortality function
"""
function HeligmanPollard(param_hp::HPParam, age::Real)
    @unpack A, B, C, D, E, F, G = param_hp
	hh = param_hp.H
    qqq = A.^((age.+B).^C) .+ D.*exp.(-E.*(log.(age) .- log.(F)).^2) .+ G.*(hh.^age)
	μ = qqq./(1 .- qqq)
	μ = max.(min.(μ, 1.0), 0.0)

    return μ
end


#plot(siler.([param], 0:110, spec = :Colchero))
#plot(siler.([param], 0:110, spec = :Scott))
#plot(siler.([param], 0:110, spec = :Bergeron))

"""
Siler survival function for individual aged aa until age tt
"""
function siler_S(param::SilerParam, aa::Real, tt::Real; spec::Symbol = :Bergeron)
	@unpack b, B, c, C, d = param
	if spec == :Colchero
		S_aa = exp( - d*aa + (1/b)*(exp(-(b*aa + B)) - exp(-B)) - (1/c)*(exp(c*aa - C) - exp(-C)) )
		S_tt = exp( - d*tt + (1/b)*(exp(-(b*tt + B)) - exp(-B)) - (1/c)*(exp(c*tt - C) - exp(-C)) )
	elseif spec == :Scott
		S_aa = exp( - d*aa + (1/b)*(exp(-b*(aa + B)) - exp(-b*B)) - (1/c)*(exp(c*(aa - C)) - exp(-c*C)) )
		S_tt = exp( - d*tt + (1/b)*(exp(-b*(tt + B)) - exp(-b*B)) - (1/c)*(exp(c*(tt - C)) - exp(-c*C)) )
	elseif spec == :Bergeron
		S_aa = exp( - d*aa + (exp(-b*(aa + B)) - exp(-b*B)) - (exp(c*(aa - C)) - exp(-c*C)) )
		S_tt = exp( - d*tt + (exp(-b*(tt + B)) - exp(-b*B)) - (exp(c*(tt - C)) - exp(-c*C)) )
	end
	if S_tt == 0.
		S_at = 0.
	else
		S_at = min.(S_tt/S_aa, 1.0)
	end
	return S_at
end


"""
Computes gradient of survivor rate from age aa to age tt, wrt θ
"""
function siler_Sgrad(param::SilerParam, aa::Real, tt::Real, θ::Symbol; spec::Symbol = :Bergeron)
	# First deriviate of surviving until aa wrt opts.param
	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		siler_S(temp_param, aa, tt, spec = spec) end, getfield(param, θ))
	return grad
end


"""
Computes remaining life expectancy from age aa given siler parameters
"""
function LE(param::SilerParam, aa::Real; spec::Symbol = :Bergeron)
    # Compute LE variable
	v = 0.
	try
		v, err = quadgk(tt -> siler_S(param, aa, tt, spec = spec),aa,505)
	catch
		v = sum(siler_S.([param], [aa], aa:505, spec = spec))
	end
    life_exp = v
    return life_exp
end

"""
Computes gradient remaining life expectancy from age aa, wrt θ
"""
function LEgrad(param::SilerParam, aa::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		LE(temp_param, aa, spec = spec) end, getfield(param, θ))
	return grad
end


"""
Computes cross-gradient remaining life expectancy from age aa, wrt θ and γ
"""
function LEcross(param::SilerParam, aa::Real, θ::Symbol, γ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, γ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, γ, x);
		LEgrad(temp_param, aa, θ, spec = spec) end, getfield(param, γ))
	return grad
end
#LEcross.([param], 0:110, [:c], [:C], spec = :Bergeron)


"""
Computes remaining lifespan inequality from age aa given siler parameters
"""
function H(param::SilerParam, aa::Real; spec::Symbol = :Bergeron)
	LE_aa = LE(param, aa, spec = spec)
	S_at = siler_S.([param], [aa], aa:500, spec = spec)
	S_at = S_at[S_at.> 0]
	H_aa = -sum(S_at.*log.(S_at))/LE_aa

	return H_aa
end


"""
Computes gradient of remaining lifespan inequality from age aa wrt θ
"""
function Hgrad(param::SilerParam, aa::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ)/0.5)(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		H(temp_param, aa; spec = spec) end, getfield(param, θ))
	return grad
end


"""
Computes remaining lifespan equality from age aa given siler parameters
"""
function h(param::SilerParam, aa::Real; spec::Symbol = :Bergeron)
	H_aa = H(param, aa, spec = spec)
	h_aa = -log(H_aa)

	return h_aa
end


"""
Computes gradient of remaining lifespan inequality from age aa wrt θ
"""
function hgrad(param::SilerParam, aa::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		h(temp_param, aa; spec = spec) end, getfield(param, θ))
	return grad
end


"""
Computes cross-gradient remaining life expectancy from age aa, wrt θ and γ
"""
function hcross(param::SilerParam, aa::Real, θ::Symbol, γ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, γ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, γ, x);
		hgrad(temp_param, aa, θ, spec = spec) end, getfield(param, γ))
	return grad
end
#hcross.([param], 0:110, [:c], [:C], spec = :Bergeron)



"""
Compute lifespan - age at a given survival probability (default 0.01)
"""
function lifespan(param; Sstar::Float64 = 0.001, spec::Symbol = :Bergeron)
	# S_t = siler_S.([param], [0], 0:200, spec = spec)
	Lst = 100.
	try
		Lst = find_zero(function(x) Sstar -  siler_S(param, 0, x, spec = spec) end,100.0)
	catch
		S_t = siler_S.([param], [0], 0:500, spec = spec)
		Lst_init = minimum(Int.(0:500)[S_t .< Sstar])
		try
			Lst = find_zero(function(x) Sstar -  siler_S(param, 0, x, spec = spec) end,Lst_init)
		catch
			display(Lst_init)
			Lst = find_zero(function(x) Sstar -  siler_S(param, 0, x, spec = spec) end,(0,  500))
		end

	end

	return Lst
end


"""
Computes gradient of remaining lifespan inequality from age aa wrt θ
"""
function lifespangrad(param::SilerParam, Sstar::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		lifespan(temp_param, Sstar = Sstar, spec = spec) end, getfield(param, θ))
	return grad
end



"""
Compute median remaining lifespan - age at with survival probability of 0.5
"""
function Lmed(param, aa; spec::Symbol = :Bergeron)
	# S_t = siler_S.([param], [aa], 0:200, spec = spec)
	Lm = 100.
	try
		Lm = find_zero(function(x) 0.5 -  siler_S(param, aa, x, spec = spec) end,100.0)
	catch
		try
			S_t = siler_S.([param], [aa], 0:500, spec = spec)
			Lm_init = minimum(Int.(0:500)[S_t .< 0.5])
			Lm = find_zero(function(x) 0.5 -  siler_S(param, aa, x, spec = spec) end,Lm_init)
		catch
			Lm = find_zero(function(x) 0.5 -  siler_S(param, aa, x, spec = spec) end,(0,500))
		end
	end
	return Lm - aa
end


"""
Computes gradient of lifespan inequality from age aa wrt θ
"""
function Lmedgrad(param::SilerParam, aa::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		Lmed(temp_param, aa, spec = spec) end, getfield(param, θ))
	return grad
end



"""
Function to calculate the ratio between mortality rates at two different ages
"""
function rμ(param::SilerParam, a1,a2; spec::Symbol = :Bergeron)
	r_out = siler(param, a1, spec = spec)/siler(param, a2, spec = spec)

	return r_out
end


"""
Compute the gradient of the ratio between mortality rates at two different ages (a1/a2)
"""
function rμgrad(param::SilerParam, a1::Real, a2::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		rμ(temp_param, a1,a2, spec = spec) end, getfield(param, θ))
	return grad
end


"""
Function to calculate the difference between mortality rates at two different ages
"""
function dμ(param::SilerParam, a1,a2; spec::Symbol = :Bergeron)
	d_out = siler(param, a1, spec = spec) - siler(param, a2, spec = spec)

	return d_out
end


"""
Compute the gradient of the difference between mortality rates at two different ages (a1/a2)
"""
function dμgrad(param::SilerParam, a1::Real, a2::Real, θ::Symbol; spec::Symbol = :Bergeron)

	grad = central_fdm(5,1,
		max_range = getfield(param, θ))(function(x)
		temp_param=deepcopy(param);
		setfield!(temp_param, θ, x);
		dμ(temp_param, a1,a2, spec = spec) end, getfield(param, θ))
	return grad
end





"""
Function to initialise an illustrative dataframe
"""
function init_illus(param::SilerParam, spec::Symbol; ages = Int.(0:110))

    illus_df = DataFrame(age = ages)
    illus_df.μ_base = siler.([param], illus_df.age, spec = spec)
    illus_df.S_base = siler_S.([param], [0.0], illus_df.age, spec = spec)
    illus_df.LE_base = LE.([param], illus_df.age, spec = spec)
    illus_df.H_base = H.([param], illus_df.age, spec = spec)
	illus_df.h_base = h.([param], illus_df.age, spec = spec)
    return illus_df
end
