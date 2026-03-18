function build_covariance_matrix(Sigma_index::Int, Sigma_param::Float64, p::Int)
    if Sigma_index == 1  # Spike
        return iszero(Sigma_param) ? Diagonal(ones(p)) : Diagonal(vcat(fill(Sigma_param, 1), ones(p - 1)))
    elseif Sigma_index == 2  # Bidiagonal
        if iszero(Sigma_param)
            return Diagonal(ones(p))
        else
            Sigma = Matrix(Bidiagonal(fill(1.0, p), fill(Sigma_param, p - 1), :L))
            Sigma[1, p] = Sigma_param
            return Sigma
        end
    else
        error("Unknown Sigma_index: $Sigma_index")
    end
end

function generate_data_matrix(dist_index::Int, Sigma::AbstractMatrix{Float64}, n::Int, p::Int)
    dist_config = DISTRIBUTIONS[dist_index]
    if dist_index == 1
        scale = 1 / dist_config.params.rate
        W = rand(Exponential(scale), n, p)
    elseif dist_index == 2
        W = rand(Chisq(dist_config.params.df), n, p)
    else
        error("Unknown dist_index")
    end
    return W * Sigma
end

function normalize_rows!(X::AbstractMatrix{T}) where {T<:AbstractFloat}
    X ./= sum(X, dims=2)
    return X
end

const SIGNIFICANCE_LEVEL = 0.95

function compute_statistics_with_estimated_parameter(X::Matrix{Float64})
    n = size(X, 1)
    p = size(X, 2)

    n1 = Int(floor(n / 2))
    n2 = n - n1
    y1 = p / n1
    y1n = p / (n1 - 1)
    y2 = p / n2
    X_stat = @view X[1:n1, :]
    X_param = @view X[(n1+1):n, :]
    B = p^2 * cov(X_stat)

    pX = p .* X_param
    fourth_central_moment = mean(x -> (x - 1)^4, pX)
    third_moment = mean(x -> x^3, pX)
    lambda_emp = mean(x -> (x - 1)^2, pX)

    nu2 = lambda_emp
    alpha1 = fourth_central_moment - 3 * lambda_emp^2
    alpha2 = -4 * lambda_emp * third_moment + 4 * lambda_emp^3 + 12 * lambda_emp^2 + 4 * lambda_emp
    h1 = -2 * third_moment + 3 * lambda_emp^2 + 5 * lambda_emp + 2
    h2 = -8 * lambda_emp * third_moment + 10 * lambda_emp^3 + 22 * lambda_emp^2 + 8 * lambda_emp
    mu2 = (1 + y1) * lambda_emp^2 + 2 * (1 + y1) * lambda_emp * h1 + y1 * (alpha1 + alpha2)
    V1 = 2 * y1 * lambda_emp^2 + y1 * (alpha1 + alpha2)
    V2 = 4 * y1 * (2 + y1) * (1 + 2 * y1) * lambda_emp^4 + 4 * y1 * (1 + y1)^2 * lambda_emp^2 * (alpha1 + alpha2)
    V12 = 2 * y1 * (1 + y1) * lambda_emp * (2 * lambda_emp^2 + alpha1 + alpha2)
    mu_T_hat = lambda_emp^2 * y1n + mu2 / p - (2 * lambda_emp * h1 + lambda_emp^2) / (p - 1)
    sigma2_T_hat = 4 * lambda_emp^2 * V1 - 4 * lambda_emp * V12 + V2
    mu_lambda_hat = -2 * lambda_emp * y1 * h1
    sigma2_lambda_hat = (2 * lambda_emp * y1)^2 * y2 * (fourth_central_moment - lambda_emp^2 + h2 - 2 * lambda_emp * h1)
    
    c = p * nu2 / (p - 1)
    c_over_p = c / p
    s = 0.0
    @inbounds for j in 1:p, i in 1:p
        diff = B[i, j] + c_over_p - (i == j ? c : 0.0)
        s += diff * diff
    end
    test_stat = s / p
    test_stat_normalized = (p * (test_stat - mu_T_hat) - mu_lambda_hat) / sqrt(sigma2_T_hat + sigma2_lambda_hat)
    decision = test_stat_normalized >= quantile(Normal(0, 1), SIGNIFICANCE_LEVEL) ? 1 : 0

    return test_stat_normalized, decision
end

function compute_statistics_with_true_parameter(X::Matrix{Float64}, dist_index::Int)
    n = size(X, 1)
    p = size(X, 2)
    y = p / n
    yn = p / (n - 1)
    B = p^2 * cov(X)
    if (dist_index == 1) # Exponential 
        lambda = 1.0
        fourth_central_moment = 9.0
        third_moment = 6.0
    elseif (dist_index == 2) # Chi-squared(df)
        dist_config = DISTRIBUTIONS[dist_index]
        df = dist_config.params.df
        mu = df
        sigma2 = 2 * df
        m1 = df
        m2 = df * (df + 2)
        m3 = df * (df + 2) * (df + 4)
        m4 = df * (df + 2) * (df + 4) * (df + 6)
        fourth_central_moment = (mu^4 - 4 * mu^3 * m1 + 6 * mu^2 * m2 - 4 * mu * m3 + m4) / mu^4
        third_moment = df * (df + 2) * (df + 4) / mu^3
        lambda = sigma2 / mu^2
    else
        error("Unknown dist_index: $dist_index")
    end
    alpha1 = fourth_central_moment - 3 * lambda^2
    alpha2 = -4 * lambda * third_moment + 4 * lambda^3 + 12 * lambda^2 + 4 * lambda
    h1 = -2 * third_moment + 3 * lambda^2 + 5 * lambda + 2
    mu2 = (1 + y) * lambda^2 + 2 * (1 + y) * lambda * h1 + y * (alpha1 + alpha2)
    nu2 = lambda + h1 / p
    V1 = 2 * y * lambda^2 + y * (alpha1 + alpha2)
    V2 = 4 * y * (2 + y) * (1 + 2 * y) * lambda^4 + 4 * y * (1 + y)^2 * lambda^2 * (alpha1 + alpha2)
    V12 = 2 * y * (1 + y) * lambda * (2 * lambda^2 + alpha1 + alpha2)
    mu_T = lambda^2 * yn + mu2 / p - (2 * lambda * h1 + lambda^2) / (p - 1)
    sigma2_T = 4 * lambda^2 * V1 - 4 * lambda * V12 + V2

    c = p * nu2 / (p - 1)
    c_over_p = c / p
    s = 0.0
    @inbounds for j in 1:p, i in 1:p
        diff = B[i, j] + c_over_p - (i == j ? c : 0.0)
        s += diff * diff
    end
    test_stat = s / p
    test_stat_normalized = p * (test_stat - mu_T) / sqrt(sigma2_T)
    decision = test_stat_normalized >= quantile(Normal(0, 1), SIGNIFICANCE_LEVEL) ? 1 : 0

    return test_stat_normalized, decision
end

function single_simulation(;dist_index::Int, Sigma_index::Int, Sigma_param::Float64, p::Int, y::Float64)
    n = Int(floor(p / y))
    Sigma = build_covariance_matrix(Sigma_index, Sigma_param, p)
    W = generate_data_matrix(dist_index, Sigma, n, p)
    normalize_rows!(W)
    stat_emp, dec_emp = compute_statistics_with_estimated_parameter(W)
    stat_theo, dec_theo = compute_statistics_with_true_parameter(W, dist_index)
    return (stat_emp=stat_emp, decision_emp=dec_emp,
            stat_theo=stat_theo, decision_theo=dec_theo)
end

function merge_all_results(; 
    dist_index_seq=[1, 2], 
    Sigma_index_seq=[1, 2], 
    y_seq=[0.5, 1.0, 1.5], 
    data_root="./data/", 
    result_types=["emp", "theo"])

    mkpath(data_root)
    for result_type in result_types
        for Sigma_index in Sigma_index_seq
            combined_df = DataFrame()
            for dist_index in dist_index_seq
                for y in y_seq
                    y_dir = joinpath(data_root, "y-" * string(y))
                    file = joinpath(y_dir, "$(result_type)_$(DISTRIBUTIONS[dist_index].name)_$(COVARIANCE[Sigma_index].name).csv")
                    try
                        df = CSV.read(file, DataFrame, header=true)
                        temp_df = copy(df)
                        temp_df.Dist = fill(DISTRIBUTIONS[dist_index].name, nrow(df))
                        temp_df.Y = fill(y, nrow(df))
                        if isempty(combined_df)
                            combined_df = temp_df
                        else
                            combined_df = vcat(combined_df, temp_df)
                        end
                    catch e
                        @warn "Cannot read file $file: $e"
                    end
                end
            end
            if nrow(combined_df) > 0
                output_file = joinpath(data_root, "$(result_type)_$(COVARIANCE[Sigma_index].name).csv")
                CSV.write(output_file, combined_df, force=true)
                println("Results merged => '$output_file'")
            else
                @warn "No data available for $(result_type) results with covariance type $(COVARIANCE[Sigma_index].name)"
            end
        end
    end
end

function merge_emp_theo(; data_root="./data/",
    Sigma_index_seq=[1, 2])
    mkpath(data_root)
    covariance_types = [COVARIANCE[Sigma_index].name for Sigma_index in Sigma_index_seq]
    for cov_type in covariance_types
        emp_file = joinpath(data_root, "emp_$cov_type.csv")
        theo_file = joinpath(data_root, "theo_$cov_type.csv")
        emp_df = CSV.read(emp_file, DataFrame)
        theo_df = CSV.read(theo_file, DataFrame)
        emp_df[!, 1:2] = Int.(emp_df[!, 1:2])
        for col in names(emp_df)
            emp_df[!, col] = string.(emp_df[!, col])
        end
        for col in names(theo_df)
            theo_df[!, col] = string.("(", theo_df[!, col], ")")
        end
        theo_df[!, 1:2] .= ""
        merged_df = DataFrame([name => String[] for name in names(emp_df)])
        for i in 1:nrow(emp_df)
            push!(merged_df, Dict(names(emp_df) .=> values(emp_df[i, :])))
            push!(merged_df, Dict(names(theo_df) .=> values(theo_df[i, :])))
        end
        pretty_table(merged_df)
        output_file = joinpath(data_root, "emp_theo_$cov_type.xlsx")
        XLSX.writetable(output_file, collect(DataFrames.eachcol(merged_df)), DataFrames.names(merged_df), overwrite=true)
        println(" '$cov_type' results merged => '$output_file'")
    end
end

