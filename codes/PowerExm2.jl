using Random
using LinearAlgebra
using Statistics
using Distributions
using Printf

const P = 300
const N_SEQ = [200, 300, 600]
const GAMMA_SEQ = [-0.20, -0.25, -0.30, -0.35]
const ALPHA_LEVEL = 0.05
const BASIS_DISTRIBUTIONS = (
    (name = "Exp(5)", dist = Exponential(1.0 / 5.0), dirichlet_shape = 1.0, lambda = 1.0, alpha_sum = 2.0),
    (name = "ChiSq(1)", dist = Chisq(1), dirichlet_shape = 0.5, lambda = 2.0, alpha_sum = 16.0),
)

function cyclic_sigma_product(Z::Matrix{Float64}, gamma::Float64)
    p = size(Z, 2)
    W = similar(Z)
    @inbounds for j in 1:(p - 1)
        @views W[:, j] .= Z[:, j] .+ gamma .* Z[:, j + 1]
    end
    @views W[:, p] .= Z[:, p] .+ gamma .* Z[:, 1]
    return W
end

# If W_j are iid Gamma(shape = a, scale = theta), then
# W_j / sum_k W_k has Dirichlet(a, ..., a) distribution.
# Thus nu2 = E[(p*x_1 - 1)^2] = (p - 1)/(p*a + 1), and
# B = p*nu2/(p - 1)*G_p = p/(p*a + 1)*G_p.
# The scale/rate cancels after row normalization.
function target_B(p::Int, dirichlet_shape::Float64)
    Gp = Matrix{Float64}(I, p, p)
    Gp .-= 1.0 / p
    return (p / (p * dirichlet_shape + 1.0)) .* Gp
end

function q_stat(W::Matrix{Float64}, B::Matrix{Float64})
    p = size(W, 2)
    X = W ./ sum(W, dims = 2)
    BpN = p^2 .* cov(X)
    return sum(abs2, BpN .- B)
end

function sigma_T(p::Int, n::Int, lambda::Float64, alpha_sum::Float64)
    c = p / n
    V1 = 2.0 * c * lambda^2 + c * alpha_sum
    V2 = 4.0 * c * (2.0 + c) * (1.0 + 2.0 * c) * lambda^4 +
         4.0 * c * (1.0 + c)^2 * lambda^2 * alpha_sum
    V12 = 2.0 * lambda * c * (1.0 + c) * (2.0 * lambda^2 + alpha_sum)
    return sqrt(4.0 * lambda^2 * V1 - 4.0 * lambda * V12 + V2)
end

function monte_carlo_difference(; p::Int, n::Int, gamma::Float64, reps::Int, basis, rng::AbstractRNG)
    B = target_B(p, basis.dirichlet_shape)
    Q0 = Vector{Float64}(undef, reps)
    Q1 = Vector{Float64}(undef, reps)

    for r in 1:reps
        Z = rand(rng, basis.dist, n, p)
        W = cyclic_sigma_product(Z, gamma)
        Q0[r] = q_stat(Z, B)
        Q1[r] = q_stat(W, B)
    end

    delta_hat = mean(Q1) - mean(Q0)
    var_hat = var(Q1; corrected = true)
    z_alpha = quantile(Normal(), 1.0 - ALPHA_LEVEL)
    rhs_hat = z_alpha * sigma_T(p, n, basis.lambda, basis.alpha_sum) +
              sqrt(ALPHA_LEVEL / (1.0 - ALPHA_LEVEL) * var_hat)
    return delta_hat - rhs_hat
end

function write_tex_row(io::IO, label::String, n::Int, row)
    @printf(io, "\t\t\t%s&(%d,%d)  & %.2f & %.2f & %.2f & %.2f  \\\\\n",
            label, P, n, row[1], row[2], row[3], row[4])
end

function write_tex_table(path::String, results)
    open(path, "w") do io
        println(io, raw"\begin{table}[htbp]")
        println(io, raw"		{\color{red}\begin{tabular}{cccccc}")
        println(io, raw"			\toprule")
        println(io, raw"			&$(p,n)\;\setminus\; \gamma$ & $-0.2$ & $-0.25$ & $-0.3$ & $-0.35$ \\\\")
        println(io, raw"			\midrule")

        exp_result = results["Exp(5)"]
        chisq_result = results["ChiSq(1)"]
        write_tex_row(io, "Exp(5)", N_SEQ[1], exp_result[1, :])
        for i in 2:length(N_SEQ)
            write_tex_row(io, "", N_SEQ[i], exp_result[i, :])
        end
        write_tex_row(io, raw"$\chi^2(1)$", N_SEQ[1], chisq_result[1, :])
        for i in 2:length(N_SEQ)
            write_tex_row(io, "", N_SEQ[i], chisq_result[i, :])
        end

        println(io, raw"			\bottomrule")
        println(io, raw"		\end{tabular}}")
        println(io, raw"		\caption{{\color{red}Difference between the two sides of \eqref{eq:Frobenius-norm-ineq} for each combination of $(p,n)$ and $\gamma$.}}")
        println(io, raw"		\label{tab:diff_covtest_Sigma2}")
        println(io, raw"	\end{table}")
    end
end

function main(; reps::Int = 2000, seed::Int = 2026)
    rng = MersenneTwister(seed)
    results = Dict{String, Matrix{Float64}}()

    for basis in BASIS_DISTRIBUTIONS
        result = Matrix{Float64}(undef, length(N_SEQ), length(GAMMA_SEQ))
        println("basis data: $(basis.name)")

        for (i, n) in enumerate(N_SEQ)
            for (j, gamma) in enumerate(GAMMA_SEQ)
                diff = monte_carlo_difference(; p = P, n = n, gamma = gamma, reps = reps, basis = basis, rng = rng)
                result[i, j] = diff
                @printf("finished (p,n) = (%d,%d), gamma = %.2f, diff = %.6f\n",
                        P, n, gamma, diff)
            end
        end

        results[basis.name] = result
    end

    output_file = joinpath(pwd(), "diff_covtest_Sigma2.tex")
    write_tex_table(output_file, results)
    println("TeX table saved to $(output_file)")
    return results
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
