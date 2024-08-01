using Turing, StatsPlots, Random, Optim, StatsBase, LinearAlgebra, Optim
using TruncatedDistributions, PDMats, Parameters, FiniteDifferences, QuadGK, Roots
using CSV, DataFrames, TableView, StatsPlots, LaTeXStrings, JLD2, ProgressMeter



"""
Parameters of Siler model
"""
@with_kw mutable struct SilerParam
    b::Float64 = 1.0 # Infant speed
    B::Float64 = 1.5 # Infant scale
    c::Float64 = 0.1 # Elderly speed
    C::Float64 = 10.0 # Elderly scale
    d::Float64 = 0.005 # constant
    σ::Float64 = 0.0 # variance
end

@with_kw mutable struct HPParam
    A::Float64 = 0.1 # Infant
    B::Float64 = 0.1 # Infant
    C::Float64 = 0.1 # Infant
    D::Float64 = 0.1 # Middle age
    E::Float64 = 10. # Middle age
    F::Float64 = 25.0 # Middle age
    G::Float64 = 0.1 # Elderly
    H::Float64 = 1.0 # Elderly
    σ::Float64 = 0.0 # Variance
end

# Various Siler functions and their derivatives
include("siler_fns.jl")
# Helper functions
include("helpers.jl")
# Turing models that estimate static and dynamic Siler models
include("models.jl")
# Functions to plot results
include("plotting.jl")
