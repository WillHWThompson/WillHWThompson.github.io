include("simple_corrected_eq.jl")

println("Generating basic data for interactive plots...")

# Nonlinearity parameters
kappa1_sign = sqrt(2/π)
kappa_star_sign = sqrt(1 - 2/π)

kappa1_erf = 2/sqrt(3*π) 
kappa_star_erf = 0.2003

# Function to safely compute generalization error
function safe_compute_error(alpha, gamma, kappa1, kappa_star, lambda_val; loss_type=:classification, Delta=0.0)
    try
        V, Q, M = corrected_fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda_val; 
                                                 loss_type=loss_type, Delta=Delta)
        gen_error = 1.0 + Q - 2*M
        
        if isnan(gen_error) || gen_error < 0 || gen_error > 10
            return NaN
        end
        return gen_error
    catch e
        return NaN
    end
end

# Generate 1D data for different scenarios
println("Generating 1D data...")

# Parameters
p_n_values = collect(0.3:0.1:2.5)
lambda_values = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
Delta_values = [0.0, 0.1, 0.5, 1.0]

# Generate and save data in simple CSV-like format
function write_data_file(filename, headers, data)
    open(filename, "w") do f
        println(f, join(headers, ","))
        for row in data
            println(f, join(row, ","))
        end
    end
end

# 1D Classification data (sign nonlinearity)
println("  Classification with sign nonlinearity...")
class_data = []
push!(class_data, ["nonlinearity", "lambda", "p_n", "n_d", "gen_error"])

for lambda_val in lambda_values
    for p_n in p_n_values
        alpha = 1.0 / p_n
        gamma = alpha / 2.0  # Fixed n/d = 2
        
        gen_error = safe_compute_error(alpha, gamma, kappa1_sign, kappa_star_sign, lambda_val)
        if !isnan(gen_error)
            push!(class_data, ["sign", lambda_val, p_n, 2.0, gen_error])
        end
    end
end

# Add erf nonlinearity classification data
for lambda_val in lambda_values
    for p_n in p_n_values
        alpha = 1.0 / p_n
        gamma = alpha / 2.0
        
        gen_error = safe_compute_error(alpha, gamma, kappa1_erf, kappa_star_erf, lambda_val)
        if !isnan(gen_error)
            push!(class_data, ["erf", lambda_val, p_n, 2.0, gen_error])
        end
    end
end

write_data_file("data/classification_1d.csv", 
                ["nonlinearity", "lambda", "p_n", "n_d", "gen_error"], 
                class_data[2:end])

# 1D Regression data (sign nonlinearity)
println("  Regression with sign nonlinearity...")
reg_data = []

for lambda_val in lambda_values
    for Delta_val in Delta_values
        for p_n in p_n_values
            alpha = 1.0 / p_n
            gamma = alpha / 2.0
            
            gen_error = safe_compute_error(alpha, gamma, kappa1_sign, kappa_star_sign, 
                                         lambda_val, loss_type=:regression, Delta=Delta_val)
            if !isnan(gen_error)
                push!(reg_data, ["sign", lambda_val, Delta_val, p_n, 2.0, gen_error])
            end
        end
    end
end

write_data_file("data/regression_1d.csv",
                ["nonlinearity", "lambda", "Delta", "p_n", "n_d", "gen_error"],
                reg_data)

# 3D Surface data  
println("  3D surface data...")
surface_data = []

p_n_3d = collect(0.4:0.3:2.2)
n_d_3d = collect(0.8:0.3:2.0)
lambda_3d = [1e-4, 1e-3, 1e-2, 1e-1]

for nonlin in ["sign", "erf"]
    kappa1 = nonlin == "sign" ? kappa1_sign : kappa1_erf
    kappa_star = nonlin == "sign" ? kappa_star_sign : kappa_star_erf
    
    for lambda_val in lambda_3d
        for n_d in n_d_3d
            for p_n in p_n_3d
                alpha = 1.0 / p_n
                gamma = alpha / n_d
                
                gen_error = safe_compute_error(alpha, gamma, kappa1, kappa_star, lambda_val)
                if !isnan(gen_error)
                    push!(surface_data, [nonlin, lambda_val, p_n, n_d, gen_error])
                end
            end
        end
    end
end

write_data_file("data/surface_3d.csv",
                ["nonlinearity", "lambda", "p_n", "n_d", "gen_error"],
                surface_data)

# Write configuration
open("data/config.txt", "w") do f
    println(f, "# Configuration for interactive plots")
    println(f, "lambda_values: ", join(lambda_values, ","))
    println(f, "Delta_values: ", join(Delta_values, ","))
    println(f, "p_n_range: 0.3,2.5")
    println(f, "n_d_range: 0.8,2.0")
    println(f, "nonlinearities: sign,erf")
end

println("Data generation complete!")
println("Files created:")
println("  - data/classification_1d.csv")
println("  - data/regression_1d.csv") 
println("  - data/surface_3d.csv")
println("  - data/config.txt")

# Validation - check for double descent spike
println("\nValidation:")
if !isempty(class_data)
    # Find spike at p/n ≈ 1 for low lambda
    spike_found = false
    for row in class_data[2:end]  # Skip header
        if row[1] == "sign" && row[2] == 1e-3 && abs(row[3] - 1.0) < 0.15
            println("Found spike at p/n = $(row[3]), λ = $(row[2]): gen_error = $(row[5])")
            spike_found = true
            if row[5] > 0.5
                println("✅ Double descent spike confirmed!")
            else
                println("⚠️  Spike value lower than expected")
            end
            break
        end
    end
    
    if !spike_found
        println("❌ No spike found around p/n = 1")
    end
else
    println("❌ No classification data generated!")
end