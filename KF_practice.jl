"""
Example Kalman Filter
"""

if occursin("jashwin", pwd())
    cd("C://Users/jashwin/Documents/GitHub/MortalityEstimation/")
else
    cd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
end

using LinearAlgebra, Statistics, Plots, Distributions


# set up prior objects
Σ = [0.4  0.3
     0.3  0.45]
x̂ = [0.2, -0.2]

# define G and R from the equation y = Gx + N(0, R)
G = I # this is a generic identity object that conforms to the right dimensions
R = 0.5 .* Σ

# define A and Q
A = [1.2  0
     0   -0.2]
Q = 0.3Σ

y = [2.3, -1.9]

# plotting objects
x_grid = range(-1.5, 2.9, length = 100)
y_grid = range(-3.1, 1.7, length = 100)

# generate distribution
dist = MvNormal(x̂, Σ)
two_args_to_pdf(dist) = (x, y) -> pdf(dist, [x, y]) # returns a function to be plotted

# Prior and observation
contour(x_grid, y_grid, two_args_to_pdf(dist), fill = false,
        color = :lighttest, cbar = false)
contour!(x_grid, y_grid, two_args_to_pdf(dist), fill = false, lw=1,
         color = :grays, cbar = false)
annotate!(y[1], y[2], "y", color = :black)


# define posterior objects
M = Σ * G' * inv(G * Σ * G' + R)
x̂_F = x̂ + M * (y - G * x̂)
Σ_F = Σ - M * G * Σ

# plot the new density on the old plot
newdist = MvNormal(x̂_F, Symmetric(Σ_F)) # because Σ_F
contour!(x_grid, y_grid, two_args_to_pdf(newdist), fill = false,
         color = :lighttest, cbar = false)
contour!(x_grid, y_grid, two_args_to_pdf(newdist), fill = false, levels = 7,
         color = :grays, cbar = false)
contour!(x_grid, y_grid, two_args_to_pdf(dist), fill = false, levels = 7, lw=1,
         color = :grays, cbar = false)


# get the predictive distribution
new_x̂ = A * x̂_F
new_Σ = A * Σ_F * A' + Q
predictdist = MvNormal(new_x̂, Symmetric(new_Σ))

# plot Density 3
contour(x_grid, y_grid, two_args_to_pdf(predictdist), fill = false, lw = 1,
    color = :lighttest, cbar = false)
contour!(x_grid, y_grid, two_args_to_pdf(dist),
    color = :grays, cbar = false)
contour!(x_grid, y_grid, two_args_to_pdf(newdist), fill = false, levels = 7,
    color = :grays, cbar = false)
annotate!(y[1], y[2], "y", color = :black)
