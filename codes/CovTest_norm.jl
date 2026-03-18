using LinearAlgebra
using Distributions
using Random
using Statistics
using ProgressMeter
using PrettyTables
using Plots
using LaTeXStrings
using Printf

cd(@__DIR__)

function normalize_rows(W::Matrix{Float64})
    row_sums = sum(W, dims=2)
    return W ./ row_sums
end

function compare_test_statistics(; 
    n_seq=[100, 200, 400, 600, 800, 1000], 
    p=200, 
    loop_num=1000, 
    )
    
    results_stat0 = zeros(loop_num, length(n_seq))
    results_stat1 = zeros(loop_num, length(n_seq))
    results_stat2 = zeros(loop_num, length(n_seq))
    Bidiagonal_value1 = 0.5
    Bidiagonal_value2 = -0.15
    for (n_idx, n) in enumerate(n_seq)
        println("n = $n, p = $p")
        
        @showprogress "Running simulations..." for t in 1:loop_num
            W0 = rand(Exponential(1), n, p)
            X0 = normalize_rows(W0)
            
            M = Bidiagonal(fill(1.0, p), fill(Bidiagonal_value1, p - 1), :L)
            Sigma1 = Matrix(M)
            Sigma1[1,p] = Bidiagonal_value1
            W1 = W0 * Sigma1
            X1 = normalize_rows(W1)

            M = Bidiagonal(fill(1.0, p), fill(Bidiagonal_value2, p - 1), :L)
            Sigma2 = Matrix(M)
            Sigma2[1,p] = Bidiagonal_value2
            W2 = W0 * Sigma2
            X2 = normalize_rows(W2)

            B0 = p^2 * cov(X0)
            B1 = p^2 * cov(X1)
            B2 = p^2 * cov(X2)
            G = Matrix{Float64}(I, p, p) - (1 / p) * ones(p, p)

            results_stat0[t, n_idx] = norm(B0 - (p / (p + 1)) * G, 2)^2 / p
            results_stat1[t, n_idx] = norm(B1 - (p / (p + 1)) * G, 2)^2 / p
            results_stat2[t, n_idx] = norm(B2 - (p / (p + 1)) * G, 2)^2 / p
        end
    end

    means_stat0 = vec(mean(results_stat0, dims=1))
    means_stat1 = vec(mean(results_stat1, dims=1))
    means_stat2 = vec(mean(results_stat2, dims=1))
    
    fig = plot(n_seq, means_stat0, 
        label=L"\frac{1}{p}\|\mathbf{B}_0 - \mathbf{B}\|_F^2", 
        marker=:circle, linewidth=2, markersize=5, color=:black, 
        xticks=(n_seq[1:2:end], string.(n_seq[1:2:end])),
        yticks=([0.5, 1, 1.5, 2], ["0.5", "1", "1.5", "2"]),
        xlabel="Sample size n",
        ylabel="Distance",
        legend=:topright,
        legendfontsize=11,
        guidefontsize=12,
        tickfontsize=10,
        fontfamily="Palatino",
        size=(700, 450),
        linestyle=:dash
    )
    plot!(n_seq, means_stat1, 
        label=L"\frac{1}{p}\|\mathbf{B}_1 - \mathbf{B}\|_F^2, (\alpha = %$(Bidiagonal_value1))",
        marker=:diamond, linewidth=2, markersize=5, color=:red)
    plot!(n_seq, means_stat2, 
        label=L"\frac{1}{p}\|\mathbf{B}_1 - \mathbf{B}\|_F^2, (\alpha = %$(Bidiagonal_value2))", 
        marker=:square, linewidth=2, markersize=5, color=:blue)

    critical_value = 7*(p+1)^2/3/p
    plot!([critical_value, critical_value], [minimum(means_stat0), 0.65], linestyle=:dash, label="", linewidth=1, color=:gray)
    annotate!([(critical_value, 
        0.8*maximum(vcat(means_stat0, means_stat1)), 
        (L"n=\frac{7(p+1)^2}{3p}\approx %$(Int(round(critical_value)))", 10, :gray))])
    annotate!(minimum(n_seq)+(critical_value-minimum(n_seq))/2, maximum(vcat(means_stat0, means_stat1))*0.3, 
        text("Detection fails\n " * L"{\scriptstyle({\color{red}{\|\mathbf{B}_1 - \mathbf{B}\|_F^2}} < {\color{black}{\|\mathbf{B}_0 - \mathbf{B}\|_F^2}})}", :black, 11, "Palatino"))
    annotate!((critical_value+maximum(n_seq))/2, maximum(vcat(means_stat0, means_stat1))*0.5, 
        text("Detection succeeds\n " * L"{\scriptstyle({\color{red}{\|\mathbf{B}_1 - \mathbf{B}\|_F^2}} > {\color{black}{\|\mathbf{B}_0 - \mathbf{B}\|_F^2}})}", :black, 11, "Palatino"))
    annotate!((critical_value+maximum(n_seq))/2, maximum(vcat(means_stat0, means_stat2))*0.5, 
        text("Detection succeeds\n " * L"{\scriptstyle({\color{blue}{\|\mathbf{B}_1 - \mathbf{B}\|_F^2}} > {\color{black}{\|\mathbf{B}_0 - \mathbf{B}\|_F^2}})}", :black, 11, "Palatino"))
    return (results_stat0=results_stat0, results_stat1=results_stat1, results_stat2=results_stat2, fig = fig)
end

Random.seed!(1234)
n_seq =200:50:1200;
p = 200;
result = compare_test_statistics(n_seq=n_seq, p=p, loop_num=1000);
display(result.fig)
savefig(result.fig, "../figure/CovTest/Frobenius_norm_comparison_nondiag.pdf")



function compare_test_statistics_diag(; 
    n_seq=[100, 200, 400, 600, 800, 1000], 
    p=200, 
    loop_num=1000
    )
    
    Sigma = Diagonal(vcat(fill(3), ones(p-1)))

    results_stat0 = zeros(loop_num, length(n_seq))
    results_stat1 = zeros(loop_num, length(n_seq))

    for (n_idx, n) in enumerate(n_seq)
        println("n = $n, p = $p")
        
        @showprogress "Running simulations..." for t in 1:loop_num
            W0 = rand(Uniform(0, 1), n, p)
            X0 = normalize_rows(W0)
            
            W1 = W0 * Sigma
            X1 = normalize_rows(W1)

            B0 = p^2 * cov(X0)
            B1 = p^2 * cov(X1)
            G = Matrix{Float64}(I, p, p) - (1 / p) * ones(p, p)

            results_stat0[t, n_idx] = norm(B0 - (p / (p + 1)) * G, 2)^2 / p
            results_stat1[t, n_idx] = norm(B1 - (p / (p + 1)) * G, 2)^2 / p
        end
    end

    means_stat0 = vec(mean(results_stat0, dims=1))
    means_stat1 = vec(mean(results_stat1, dims=1))
    
    fig = plot(n_seq, means_stat0, 
        label=L"\frac{1}{p}\|\mathbf{B}_0 - \mathbf{B}\|_F^2", 
        marker=:circle, linewidth=2, markersize=5, 
        xticks=(n_seq[1:2:end], string.(n_seq[1:2:end])),
        xlabel="Sample size n",
        ylabel="Distance",
        legend=:topright,
        legendfontsize=11,
        guidefontsize=12,
        tickfontsize=10,
        fontfamily="Palatino",
        size=(700, 450)
    )
    plot!(n_seq, means_stat1, 
        label=L"\frac{1}{p}\|\mathbf{B}_1 - \mathbf{B}\|_F^2", 
        marker=:square, linewidth=2, markersize=5)
    annotate!((700, 0.53, 
        text("Detection succeeds\n " * L"{({\color{red}{\|\mathbf{B}_1 - \mathbf{B}\|_F^2}} > {\color{blue}{\|\mathbf{B}_0 - \mathbf{B}\|_F^2}})}", :black, 12, "Palatino")))

    return (results_stat0=results_stat0, results_stat1=results_stat1, fig = fig)
end


Random.seed!(1234)
n_seq =200:50:1200;
p = 200;
result = compare_test_statistics_diag(n_seq=n_seq, p=p, loop_num=1000);
display(result.fig)
savefig(result.fig, "../figure/CovTest/Frobenius_norm_comparison_diag.pdf")