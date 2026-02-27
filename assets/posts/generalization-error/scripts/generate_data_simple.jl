using JSON

# Import our corrected fixed point iteration functions  
include("simple_corrected_eq.jl")

println("Generating data for interactive plots...")

# Define nonlinearity parameters
nonlinearities = Dict(
    "sign" => (kappa1=sqrt(2/π), kappa_star=sqrt(1-2/π)),
    "erf" => (kappa1=2/sqrt(3*π), kappa_star=0.2003),
    "identity" => (kappa1=1.0, kappa_star=0.0)
)

# Function to safely compute generalization error
function safe_compute_error(alpha, gamma, kappa1, kappa_star, lambda_val; loss_type=:classification, Delta=0.0)
    try
        V, Q, M = corrected_fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda_val; 
                                                 loss_type=loss_type, Delta=Delta)
        gen_error = 1.0 + Q - 2*M
        
        if isnan(gen_error) || gen_error < 0 || gen_error > 10
            return nothing
        end
        return gen_error
    catch e
        return nothing
    end
end

# Generate 1D data (p/n curves)
println("Generating 1D curves...")
p_n_values = collect(0.2:0.1:3.0)
lambda_values = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
Delta_values = [0.0, 0.1, 0.5, 1.0]

data_1d = Dict()

for nonlin_name in keys(nonlinearities)
    params = nonlinearities[nonlin_name] 
    data_1d[nonlin_name] = Dict()
    println("  Processing $nonlin_name...")
    
    for lambda_val in lambda_values
        key_lambda = "lambda_$(lambda_val)"
        data_1d[nonlin_name][key_lambda] = Dict()
        
        # Classification data
        p_n_array = Float64[]
        gen_error_array = Float64[]
        
        for p_n in p_n_values
            alpha = 1.0 / p_n
            gamma = alpha / 2.0  # Fixed n/d = 2
            
            gen_error = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star, lambda_val)
            if gen_error !== nothing
                push!(p_n_array, p_n)
                push!(gen_error_array, gen_error)
            end
        end
        
        data_1d[nonlin_name][key_lambda]["classification"] = Dict(
            "p_n" => p_n_array,
            "gen_error" => gen_error_array
        )
        
        # Regression data for different noise levels
        data_1d[nonlin_name][key_lambda]["regression"] = Dict()
        
        for Delta_val in Delta_values
            key_delta = "Delta_$(Delta_val)"
            p_n_reg = Float64[]
            gen_error_reg = Float64[]
            
            for p_n in p_n_values
                alpha = 1.0 / p_n
                gamma = alpha / 2.0
                
                gen_error = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star, 
                                             lambda_val, loss_type=:regression, Delta=Delta_val)
                if gen_error !== nothing
                    push!(p_n_reg, p_n)
                    push!(gen_error_reg, gen_error)
                end
            end
            
            data_1d[nonlin_name][key_lambda]["regression"][key_delta] = Dict(
                "p_n" => p_n_reg,
                "gen_error" => gen_error_reg
            )
        end
    end
end

# Generate 3D surface data
println("Generating 3D surface data...")
p_n_3d = collect(0.3:0.2:2.5)
n_d_3d = collect(0.7:0.2:2.3)

data_3d = Dict()

for nonlin_name in ["sign", "erf"]  # Main nonlinearities for 3D
    params = nonlinearities[nonlin_name]
    data_3d[nonlin_name] = Dict()
    println("  Processing $nonlin_name for 3D...")
    
    for lambda_val in [1e-4, 1e-3, 1e-2, 1e-1]
        key_lambda = "lambda_$(lambda_val)"
        
        p_n_grid = Float64[]
        n_d_grid = Float64[]  
        gen_error_grid = Float64[]
        
        for n_d in n_d_3d
            for p_n in p_n_3d
                alpha = 1.0 / p_n
                gamma = alpha / n_d
                
                gen_error = safe_compute_error(alpha, gamma, params.kappa1, params.kappa_star, lambda_val)
                if gen_error !== nothing
                    push!(p_n_grid, p_n)
                    push!(n_d_grid, n_d)
                    push!(gen_error_grid, gen_error)
                end
            end
        end
        
        data_3d[nonlin_name][key_lambda] = Dict(
            "p_n" => p_n_grid,
            "n_d" => n_d_grid,
            "gen_error" => gen_error_grid
        )
    end
end

# Save data as JSON
println("Saving data...")

open("data/generalization_data_1d.json", "w") do f
    JSON.print(f, data_1d, 2)
end

open("data/generalization_data_3d.json", "w") do f
    JSON.print(f, data_3d, 2)
end

# Create configuration file
config = Dict(
    "nonlinearities" => collect(keys(nonlinearities)),
    "lambda_values" => lambda_values,
    "Delta_values" => Delta_values,
    "p_n_range" => [0.2, 3.0],
    "n_d_range" => [0.7, 2.3],
    "description" => "Interactive generalization error data"
)

open("data/config.json", "w") do f
    JSON.print(f, config, 2)
end

println("Data generation complete!")
println("Files created:")
println("  - data/generalization_data_1d.json")
println("  - data/generalization_data_3d.json") 
println("  - data/config.json")

# Validation
println("\nValidation:")
sign_data = data_1d["sign"]["lambda_0.001"]["classification"]
if !isempty(sign_data["p_n"])
    # Find point closest to p/n = 1
    idx = argmin(abs.(sign_data["p_n"] .- 1.0))
    spike_value = sign_data["gen_error"][idx]
    p_n_spike = sign_data["p_n"][idx]
    println("Spike at p/n ≈ $p_n_spike: gen_error = $spike_value")
    println("✅ Expected spike > 0.5 for double descent" * (spike_value > 0.5 ? " - FOUND!" : " - NOT FOUND"))
else
    println("❌ No data generated!")
end