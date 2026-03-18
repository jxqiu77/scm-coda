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

include("CoDA_utils.jl")

const DISTRIBUTIONS = Dict(
    1 => (name="Exp", params=(rate=5,)), # Exponential
    2 => (name="ChiSq", params=(df=1,)) # Chi-squared
)

const COVARIANCE = Dict(
    1 => (name="Spike", params=([1.0, 3.0, 3.5, 4.0, 4.5])),
    2 => (name="Bidiagonal", params=([0.0, -0.20, -0.25, -0.30, -0.35]))
)

function run_batch_simulation(; dist_index_seq::Vector{Int}=sort(collect(keys(DISTRIBUTIONS))),
    Sigma_index_seq::Vector{Int}=sort(collect(keys(COVARIANCE))),
    y_seq::Vector{Float64}=[0.5, 1.0, 1.5],
    dim_seq::Vector{Int}=[150, 300, 450, 600],
    loop_num::Int=2000,
    data_root::String="../data/CovTest/")

    Sigma_param_seq_all = Dict(i => COVARIANCE[i].params for i in keys(COVARIANCE))
    data_dirs = Dict(y => joinpath(data_root, "y-" * string(y)) for y in y_seq)
    for Sigma_index in Sigma_index_seq
        Sigma_param_seq = Sigma_param_seq_all[Sigma_index]
        for dist_index in dist_index_seq
            for y in y_seq
                data_dir = data_dirs[y]
                result_theo = zeros(length(dim_seq), length(Sigma_param_seq) + 2)
                result_theo[:, 1] = dim_seq
                result_theo[:, 2] = Int.(dim_seq ./ y)
                result_emp = copy(result_theo)
                for (dim_index, p) in enumerate(dim_seq)
                    n = Int(p / y)
                    for (Sigma_param_index, Sigma_param) in enumerate(Sigma_param_seq)
                        decisions_theo = zeros(loop_num)
                        decisions_emp = zeros(loop_num)
                        Random.seed!(1234) 
                        progress = Progress(
                            loop_num; 
                            desc = "Simulating $(loop_num) runs",
                            barglyphs = BarGlyphs('=', '=', '>', '.', ']')
                        )
                        for t in 1:loop_num
                            result = single_simulation(; dist_index = dist_index, Sigma_index = Sigma_index, Sigma_param = Sigma_param, p = p, y = y)
                            decisions_emp[t] = result.decision_emp
                            decisions_theo[t] = result.decision_theo 
                            next!(progress)
                        end
                        result_emp[dim_index, Sigma_param_index+2] = mean(decisions_emp)
                        result_theo[dim_index, Sigma_param_index+2] = mean(decisions_theo)
                        println("dist: $(DISTRIBUTIONS[dist_index].name), Sigma: $(COVARIANCE[Sigma_index].name), (p, n) = ($p, $(Int(p/y)))")
                        println("Empirical:")
                        pretty_table(result_emp; column_labels=[:p, :n, Symbol.(Sigma_param_seq)...])
                        println("Theoretical:")
                        pretty_table(result_theo; column_labels=[:p, :n, Symbol.(Sigma_param_seq)...])
                    end
                end
                mkpath(data_dir)
                df_emp = DataFrame(result_emp, [:p, :n, Symbol.(Sigma_param_seq)...])
                CSV.write(joinpath(data_dir, "emp_$(DISTRIBUTIONS[dist_index].name)_$(COVARIANCE[Sigma_index].name).csv"), df_emp, writeheader=true)
                df_theo = DataFrame(result_theo, [:p, :n, Symbol.(Sigma_param_seq)...])
                CSV.write(joinpath(data_dir, "theo_$(DISTRIBUTIONS[dist_index].name)_$(COVARIANCE[Sigma_index].name).csv"), df_theo, writeheader=true)
            end
        end
    end
end

cd(@__DIR__)
Sigma_index_seq = [1, 2]
loop_num = 2000
run_batch_simulation(dist_index_seq=[1, 2], Sigma_index_seq=Sigma_index_seq, y_seq=[0.5, 1.0, 1.5], dim_seq=[150, 300, 450, 600], loop_num=loop_num, data_root="../data/CovTest/")
merge_all_results(dist_index_seq=[1, 2], Sigma_index_seq=Sigma_index_seq, y_seq=[0.5, 1.0, 1.5], data_root="../data/CovTest/")
merge_emp_theo(; Sigma_index_seq=Sigma_index_seq, data_root="../data/CovTest/")
