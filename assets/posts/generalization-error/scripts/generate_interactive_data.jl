using CSV
using DataFrames
using JSON

# Import our corrected fixed point iteration functions
include("simple_corrected_eq.jl")

println("Generating comprehensive dataset for interactive plots...")

# Define parameter ranges for comprehensive coverage
lambda_values = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
p_n_values = collect(0.1:0.05:3.0)  # p/n ratios
n_d_values = collect(0.5:0.1:2.5)   # n/d ratios  
Delta_values = [0.0, 0.01, 0.1, 0.5, 1.0]  # Noise levels

# Define nonlinearity parameters
nonlinearities = Dict(
    "sign" => (kappa1=sqrt(2/π), kappa_star=sqrt(1-2/π)),
    "erf" => (kappa1=2/sqrt(3*π), kappa_star=0.2003),
    "identity" => (kappa1=1.0, kappa_star=0.0),
    "tanh" => (kappa1=0.6, kappa_star=0.4)  # Approximate values
)

# Function to safely compute generalization error with error handling
function safe_compute_error(alpha, gamma, kappa1, kappa_star, lambda_val; loss_type=:classification, Delta=0.0)
    try
        V, Q, M = corrected_fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda_val; 
                                                 loss_type=loss_type, Delta=Delta)
        rho = 1.0
        gen_error = rho + Q - 2*M
        
        # Sanity check: ensure reasonable values
        if isnan(gen_error) || gen_error < 0 || gen_error > 10
            return NaN
        end
        return gen_error
    catch e
        println("Error at α=$alpha, γ=$gamma, λ=$lambda_val: $e")
        return NaN
    end
end

# Generate 1D data (p/n vs generalization error) for different parameters
println("Generating 1D data (p/n curves)...")
data_1d = DataFrame()

for nonlin_name in keys(nonlinearities)
    params = nonlinearities[nonlin_name]
    println("  Processing $nonlin_name nonlinearity...")
    
    for lambda_val in lambda_values
        for Delta_val in Delta_values
            # Fixed n/d = 2 for 1D curves
            n_d_fixed = 2.0
            
            for (i, p_n) in enumerate(p_n_values)
                alpha = 1.0 / p_n
                gamma = alpha / n_d_fixed
                
                # Classification 
                gen_error_class = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star, 
                                                   lambda_val, loss_type=:classification, Delta=0.0)
                
                # Regression
                gen_error_reg = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star,
                                                 lambda_val, loss_type=:regression, Delta=Delta_val)
                
                push!(data_1d, (
                    nonlinearity=nonlin_name,
                    lambda=lambda_val, 
                    Delta=Delta_val,
                    p_n=p_n,
                    n_d=n_d_fixed,
                    alpha=alpha,
                    gamma=gamma,
                    gen_error_classification=gen_error_class,
                    gen_error_regression=gen_error_reg
                ))
                
                # Progress indicator
                if i % 10 == 0
                    println("    p/n = $p_n")
                end
            end
        end
    end
end

# Save 1D data
CSV.write("data/generalization_error_1d.csv", data_1d)
println("Saved 1D data to data/generalization_error_1d.csv")

# Generate 3D surface data (n/d vs p/n vs generalization error)
println("Generating 3D surface data...")
data_3d = DataFrame()

# Use fewer points for 3D to keep file size reasonable
p_n_3d = collect(0.2:0.15:2.5)  
n_d_3d = collect(0.6:0.15:2.2)

for nonlin_name in ["sign", "erf"]  # Focus on main nonlinearities for 3D
    params = nonlinearities[nonlin_name]
    println("  Processing $nonlin_name for 3D surface...")
    
    for lambda_val in [1e-4, 1e-3, 1e-2, 1e-1]  # Fewer lambda values for 3D
        for (i, n_d) in enumerate(n_d_3d)
            for (j, p_n) in enumerate(p_n_3d)
                alpha = 1.0 / p_n  
                gamma = alpha / n_d
                
                # Classification only for 3D plots
                gen_error = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star,
                                             lambda_val, loss_type=:classification, Delta=0.0)
                
                push!(data_3d, (
                    nonlinearity=nonlin_name,
                    lambda=lambda_val,
                    p_n=p_n,
                    n_d=n_d,
                    alpha=alpha, 
                    gamma=gamma,
                    gen_error=gen_error
                ))
            end
            
            if i % 3 == 0
                println("    n/d = $n_d")
            end
        end
    end
end

# Save 3D data
CSV.write("data/generalization_error_3d.csv", data_3d)
println("Saved 3D data to data/generalization_error_3d.csv")

# Create metadata file for the interactive plots
metadata = Dict(
    "nonlinearities" => collect(keys(nonlinearities)),
    "lambda_range" => [minimum(lambda_values), maximum(lambda_values)],
    "lambda_values" => lambda_values,
    "Delta_range" => [minimum(Delta_values), maximum(Delta_values)],
    "Delta_values" => Delta_values,
    "p_n_range" => [minimum(p_n_values), maximum(p_n_values)],
    "n_d_range" => [minimum(n_d_values), maximum(n_d_values)],
    "description" => "Generalization error data for interactive visualization",
    "generation_date" => string(now())
)

open("data/metadata.json", "w") do f
    JSON.print(f, metadata, 2)
end

println("Saved metadata to data/metadata.json")
println("Data generation complete!")

# Quick validation - check some key results
println("\nValidation - checking for double descent spikes:")
validation_data = filter(row -> row.nonlinearity == "sign" && row.lambda == 1e-3 && 
                        abs(row.p_n - 1.0) < 0.1, data_1d)
if !isempty(validation_data)
    spike_error = validation_data[1, :gen_error_classification]
    println("Spike at p/n ≈ 1.0: $spike_error (should be > 0.5 for low λ)")
else
    println("No validation data found!")
end