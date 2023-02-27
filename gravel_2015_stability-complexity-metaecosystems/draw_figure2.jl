# Reproduce Figure 2 of Gravel et al. 2015 (10.1038/ncomms12457).
# --------------------

using Distributions
using LaTeXStrings
using LinearAlgebra
using Plots
using Plots.PlotMeasures: mm

# We begin by defining few functions to generate the different matrices needed.

"""
    create_single_interaction_matrix(S, C, mu, sigma)

Create an interaction matrix between `S` species
with a connectance of interspecific links `C`.
Interaction intensenties are drawn in a normal law of
mean `mu` and standard deviation `sigma`.
The diagonal elements of the interaction matrix are set to 0.
"""
function create_single_interaction_matrix(S, C, mu, sigma)
    A = zeros(Float64, S, S)
    for i in 1:S, j in 1:S
        i != j && rand() < C && (A[i, j] = rand(Normal(mu, sigma)))
    end
    A
end

"""
    create_independent_interaction_matrix(n_patch, S, C, mu, sigma)

Create a meta interaction matrix, a block diagonal matrix
of size (`n_patch*S`, `n_patch*S`)
where each diagonal block of size (`S`, `S`)
is an interaction matrix of a given patch.
Here each interaction matrix are drawn independently from each other.
"""
function create_independent_interaction_matrix(n_patch, S, C, mu, sigma)
    size_tot = n_patch * S
    A_meta = zeros(Float64, size_tot, size_tot)
    for i in 1:n_patch
        A = create_single_interaction_matrix(S, C, mu, sigma)
        from = 1 + (i - 1) * S
        to = from + S - 1
        A_meta[from:to, from:to] = A
    end
    A_meta
end

"""
    create_identic_interaction_matrix(n_patch, S, C, mu, sigma)

Create a meta interaction matrix, a block diagonal matrix
of size (`n_patch*S`, `n_patch*S`)
where each diagonal block of size (`S`, `S`)
is an interaction matrix of a given patch.
Here each interaction matrix are identic.
"""
function create_identic_interaction_matrix(n_patch, S, C, mu, sigma)
    A = create_single_interaction_matrix(S, C, mu, sigma)
    block_diagonal(A, n_patch)
end


"""
    block_diagonal(block, n_block)

Create a block diagonal matrix that repeat `n_block` times the `block`
in the diagonal.
"""
function block_diagonal(block, n_block)
    block_length = size(block, 1)
    size_tot = block_length * n_block
    matrix_block_diagonal = zeros(Float64, size_tot, size_tot)
    for i in 1:n_block
        from = 1 + (i - 1) * block_length
        to = from + block_length - 1
        matrix_block_diagonal[from:to, from:to] = block
    end
    matrix_block_diagonal
end

"""
    connectivity_matrix_full(n_patch)

Create a connectivity matrix between `n_patch`
such that they are all connected with each other.
"""
function connectivity_matrix_full(n_patch)
    connectivity_matrix = zeros(Float64, n_patch, n_patch)
    for i in 1:n_patch, j in 1:n_patch
        connectivity_matrix[i, j] = i == j ? -1 : 1 / (n_patch - 1)
    end
    connectivity_matrix
end

"""
    connectivity_matrix_circle(n_patch)

Create a connectivity matrix between `n_patch`
such that patch n is connected with patch n-1 and n+1,
with periodic conditions for the first and last patch.
"""
function connectivity_matrix_circle(n_patch)
    connectivity_matrix = zeros(Float64, n_patch, n_patch) - I
    for i in 1:n_patch
        off_index = i + 1 <= n_patch ? i + 1 : 1
        connectivity_matrix[i, off_index] = 0.5
        connectivity_matrix[off_index, i] = 0.5
    end
    connectivity_matrix
end

"""
    create_dispersal_matrix(connectivity_matrix, dispersal_rates)

Given a `connectivity_matrix` of size (`n_patch`, `n_patch`) and
a vector of `dispersal_rates` for each species,
return the corresponding dispersal matrix.
"""
function create_dispersal_matrix(connectivity_matrix, dispersal_rates)
    n_patch = size(connectivity_matrix, 1)
    S = length(dispersal_rates)
    size_tot = n_patch * S
    D = zeros(Float64, size_tot, size_tot)
    for patch_i in 1:n_patch, patch_j in 1:n_patch
        connectivity = connectivity_matrix[patch_i, patch_j]
        for (sp, d) in enumerate(dispersal_rates)
            row = S * (patch_i - 1) + sp
            col = S * (patch_j - 1) + sp
            D[row, col] = connectivity * d
        end
    end
    D
end

# Set up system parameters taken from Gravel et al. 2015 (cf. Figure 2).
S = 100 # species richness
n_patch = 20
C = 0.3 # connectance
sigma = 1
m = 2 # self-regulation (identic for all species and all patches)
A_metacom = create_independent_interaction_matrix(n_patch, S, C, 0, sigma)
connectivity_matrix = connectivity_matrix_full(n_patch)
plot_list = []
for d in ([0, 1, 8])
    dispersal_rates = fill(d, S)
    D = create_dispersal_matrix(connectivity_matrix, dispersal_rates)
    J = A_metacom + D - m * I # compute community matrix

    # Plot eigen values distribution.
    eigen_values = eigvals(J)
    re_eigvals = real.(eigen_values)
    im_eigvals = imag.(eigen_values)
    temp_plot = scatter(
        re_eigvals,
        im_eigvals;
        label = false,
        color = :black,
        size = (500, 500),
        markersize = 2,
    )
    xlabel!(L"\Re")
    ylabel!(L"\Im")
    title!(L"d = %$d")
    push!(plot_list, temp_plot)
end
window_size = 300
plot(
    plot_list...;
    dpi = 300,
    layout = (1, 3),
    size = (window_size * 3, window_size),
    left_margin = 3mm,
    bottom_margin = 3mm,
)
savefig(
    "gravel_2015_stability-complexity-metaecosystems/figure2_eigenvalue-distribution.png",
)
