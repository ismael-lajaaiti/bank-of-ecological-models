# Import dependencies
using Plots, LinearAlgebra, Distributions
include("utils.jl")

N = 3000 # number of mussels
landscape_size = 50 # landscape length 
ntime = 300 # total time  

# Parameters of the function relating density to movement speed
P1 = 100.0
P2 = -80
P3 = 3
D1 = 2 # 2 - Size of the direct neighborhood
D2 = 6 # 6 - Size of the cluster neighborhood
Diagonal_mat = convert(Matrix{Float64}, I(N))

# Random initial position
X = rand(N) * landscape_size
Y = rand(N) * landscape_size

# Allocating for final position
X_n, Y_n = zeros(ntime, N), zeros(ntime, N)
X_n[1, :] = X
Y_n[1, :] = Y

# Distance between each mussel in 
distance_mussel = function (X, Y)
    Δx = X .- X'
    Δy = Y .- Y'
    return sqrt.((Δx .^ 2 .+ Δy .^ 2))
end

for t in 2:(ntime-1)
    # For each mussel, we need to check the number of mussels within the range D1 or D2
    dist_t = distance_mussel(X_n[t-1, :], Y_n[t-1, :])
    # Checking which mussels are within ranges D1 and D2
    range_D1 = convert(Matrix{Float64}, dist_t .< D1) .- Diagonal_mat
    range_D2 = convert(Matrix{Float64}, dist_t .< D2) .- Diagonal_mat
    # Calculating, for each mussel, the density of the local neighborhood (area = pi r²)
    nb_D1 = sum(range_D1; dims = 2) ./ (pi * D1^2)
    nb_D2 = sum(range_D2; dims = 2) ./ (pi * D2^2)
    # β is determined from the data
    β = 1 ./ [((max(0.02, P1 * nb_D1[x] + P2 * nb_D2[x])) + P3) for x in eachindex(1:N)]
    Stepsize = -β .* log.(rand(N)) #exponential law see the paper
    Angle = rand(N) .* 360 #random direction
    coord_x = X_n[t-1, :] .+ sin.((Angle ./ 180) .* pi) .* Stepsize
    coord_y = Y_n[t-1, :] .+ cos.((Angle ./ 180) .* pi) .* Stepsize
    X_n[t, :] = coord_x
    Y_n[t, :] = coord_y
end

anim_mussel = @animate for i in collect(1:10:size(X_n)[1])
    scatter(
        X_n[i, :],
        Y_n[i, :];
        xlabel = "x",
        ylabel = "y",
        title = "Mussels, (t = " * repr(i) * ")",
        ylim = (0, 50),
        xlim = (0, 50),
    )
end

check_dir("figures")
gif(anim_mussel, "figures/mussels_dynamics.gif"; fps = 5)
