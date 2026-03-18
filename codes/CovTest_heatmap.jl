using LinearAlgebra
using Distributions
using Random
using Statistics
using DataFrames
using CSV
using ProgressMeter
using PrettyTables
using Plots
using XLSX
using LaTeXStrings

include("CoDA_utils.jl")

const DISTRIBUTIONS = Dict(
    1 => (name="Exp", params=(rate=5,)), # Exponential
    2 => (name="ChiSq", params=(df=1,)) # Chi-squared
)

const COVARIANCE = Dict(
    1 => (name="Spike", params=([1.0, 3.0, 3.5, 4.0, 4.5])),
    2 => (name="Bidiagonal", params=([0.0, -0.1, -0.15, -0.20, -0.25, -0.30, -0.35]))
)

p_seq = 100:1:130
n_seq = 240:2:330
alpha = 0.5
dist_index = 1 # Exponential
Sigma_index = 2 # Bidiagonal
loop_num = 1000 # replications for each (p, n) pair
results_df = DataFrame(p=Int[], n=Int[], power=Float64[])

progress = Progress(
        length(p_seq) * length(n_seq);
        desc = "Simulating",
        barglyphs = BarGlyphs('=', '=', '>', '.', ']')
    )
for p in p_seq
    for n in n_seq
        decision_theo = zeros(loop_num)
        for t in 1:loop_num
            Sigma = build_covariance_matrix(Sigma_index, alpha, p)
            W = generate_data_matrix(dist_index, Sigma, n, p)
            normalize_rows!(W)
            _, dec_theo = compute_statistics_with_true_parameter(W, dist_index)
            decision_theo[t] = dec_theo
        end
        push!(results_df, (p=p, n=n, power=mean(decision_theo)))
        next!(progress)
    end
end
# Save results to CSV
csv_path = "../data/CovTest/"
mkpath(csv_path)
csv_name = joinpath(csv_path, "power_results.csv")
CSV.write(csv_name, results_df)

# Load results and create heatmap
results_from_csv = CSV.read(csv_name, DataFrame)
pivot_df = unstack(results_from_csv, :p, :n, :power)
xvals = p_seq
yvals = n_seq
zvals = Matrix(pivot_df[:, Not(:p)])' 

default(fontfamily="Palatino",
    linewidth=2,
    grid=true,
    minorticks=true,
    guidefontsize=10,
    tickfontsize=8,
    legendfontsize=8)

heatmap(
    xvals,
    yvals,
    zvals;
    xlabel="p",
    ylabel="n",
    colorbar_title="Power",
    title="Power Heatmap",
    color = :cividis,
    clims = (0, 1)
)
y_func = [(2*(alpha^2 + alpha + 1) / (3*alpha)) * (p+1)^2 / p for p in p_seq]
plot!(xvals, y_func; lw=2, c=:red, lines=:dash, label=false)
idx_arrow = floor(Int, length(xvals) * 0.7)
x_arrow = xvals[idx_arrow]
y_arrow = y_func[idx_arrow]
annotate!(x_arrow + 3, y_arrow - 30, text(L"n_{\mathrm{crit}}(p,\alpha)", 13, :red))
plot!([x_arrow + 2.7, x_arrow], [y_arrow - 27, y_arrow];
c=:red, lw=1.8, arrow=:closed, label=false)

epsilon = 0.05
plot!(xvals, y_func.*(1+epsilon); lw=2, c=:red, lines=:dot, label=false)
idx_arrow = floor(Int, length(xvals) * 0.4)
x_arrow = xvals[idx_arrow] 
y_arrow = y_func[idx_arrow] + 12
annotate!(x_arrow - 5, y_arrow + 25, text(L"(1+\varepsilon)n_{\mathrm{crit}}(p,\alpha)", 13, :red))
plot!([x_arrow - 4.7, x_arrow], [y_arrow + 20, y_arrow];
c=:red, lw=1.8, arrow=:closed, label=false)

ylims!(minimum(n_seq), maximum(n_seq))

fig_path = "../figure/CovTest/"
mkpath(fig_path)
savefig(joinpath(fig_path, "power_heatmap.pdf"))
