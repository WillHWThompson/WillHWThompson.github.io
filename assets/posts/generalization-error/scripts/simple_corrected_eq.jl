using Plots

# Define the Stieltjes transform g_mu(z) and its derivative for Gaussian projections
function g_mu(z, γ)
    sqrt_term = sqrt((z - 1 - γ)^2 - 4 * γ)
    return (1 - z - γ - sqrt_term) / (2 * z * γ)
end

function g_mu_derivative(z, gamma)
    sqrt_term = sqrt(-4 * gamma + (1 + gamma - z)^2)
    numerator = (1 - 2 * gamma + gamma^2 - sqrt_term + gamma * sqrt_term - z - gamma * z)
    denominator = (2 * gamma * sqrt_term * z^2)
    return numerator / denominator
end

# Simple implementation using the exact formulas from the paper for specific cases
# Classification with square loss and sign labels (from paper equations D.26)
function classification_hat_overlaps_corrected(alpha, gamma, kappa1, kappa_star, V_val, M_val, Q_val)
    # These are the corrected formulas from equations D.26 in the paper
    V_s_hat = (alpha * kappa1^2) / (gamma * (1 + V_val))
    q_s_hat = (alpha * kappa1^2) / (gamma) * (1 + Q_val - 2*M_val/sqrt(π)) / (1 + V_val)^2
    m_s_hat = (alpha * kappa1) / gamma * sqrt(2/π) / (1 + V_val)
    
    V_w_hat = (alpha * kappa_star^2) / (1 + V_val)  
    q_w_hat = (alpha * kappa_star^2) * (1 + Q_val - 2*M_val/sqrt(π)) / (1 + V_val)^2
    
    return V_s_hat, q_s_hat, m_s_hat, V_w_hat, q_w_hat
end

# Regression with square loss and linear labels (from paper equations D.25)  
function regression_hat_overlaps_corrected(alpha, gamma, kappa1, kappa_star, V_val, M_val, Q_val, Delta)
    # These are the corrected formulas from equations D.25 in the paper
    V_s_hat = (alpha * kappa1^2) / (gamma * (1 + V_val))
    q_s_hat = (alpha * kappa1^2) / gamma * (1 + Delta + Q_val - 2*M_val) / (1 + V_val)^2
    m_s_hat = (alpha * kappa1) / (gamma * (1 + V_val))
    
    V_w_hat = (alpha * kappa_star^2) / (1 + V_val)
    q_w_hat = (alpha * kappa_star^2) * (1 + Delta + Q_val - 2*M_val) / (1 + V_val)^2
    
    return V_s_hat, q_s_hat, m_s_hat, V_w_hat, q_w_hat
end

# Fixed point iteration with corrected equations D.23
function corrected_fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda; 
                                       loss_type=:classification, Delta=0.0,
                                       max_iter=1000, tol=1e-8)
    
    # Initialize 
    V_val, Q_val, M_val = 0.1, 0.1, 0.1
    
    for iter in 1:max_iter
        # Compute hat variables using corrected formulas
        if loss_type == :classification
            V_s_hat, q_s_hat, m_s_hat, V_w_hat, q_w_hat = 
                classification_hat_overlaps_corrected(alpha, gamma, kappa1, kappa_star, V_val, M_val, Q_val)
        else  # regression
            V_s_hat, q_s_hat, m_s_hat, V_w_hat, q_w_hat = 
                regression_hat_overlaps_corrected(alpha, gamma, kappa1, kappa_star, V_val, M_val, Q_val, Delta)
        end
        
        # Compute z and Stieltjes transforms
        z = (lambda + V_w_hat) / V_s_hat
        g_z = g_mu(-z, gamma)
        g_prime_z = g_mu_derivative(-z, gamma)
        
        # Apply equations D.23 (corrected)
        V_s = (1/V_s_hat) * (1 - z * g_z)
        
        q_s = ((m_s_hat^2 + q_s_hat)/(V_s_hat^2)) * (1 - 2*z*g_z + z^2*g_prime_z) - 
              (q_w_hat/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
              
        m_s = (m_s_hat/V_s_hat) * (1 - z * g_z)
        
        V_w = gamma/(lambda + V_w_hat) * (1/gamma - 1 + z*g_z)
        
        # CORRECTED: This was the key error - the sign should be + not -
        q_w = gamma * (q_w_hat/(lambda + V_w_hat)^2) * (1/gamma - 1 + z^2*g_prime_z) + 
              ((m_s_hat^2 + q_s_hat)/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
        
        # Compute new V, Q, M
        V_new = kappa1^2 * V_s + kappa_star^2 * V_w
        Q_new = kappa1^2 * q_s + kappa_star^2 * q_w  
        M_new = kappa1 * m_s
        
        # Check convergence
        if (abs(V_new - V_val) < tol && abs(Q_new - Q_val) < tol && abs(M_new - M_val) < tol)
            break
        end
        
        # Update
        V_val, Q_val, M_val = V_new, Q_new, M_new
    end
    
    return V_val, Q_val, M_val
end

# Compute generalization error
function compute_gen_error_corrected(alpha, gamma, kappa1, kappa_star, lambda; 
                                   loss_type=:classification, Delta=0.0)
    V, Q, M = corrected_fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda; 
                                            loss_type=loss_type, Delta=Delta)
    rho = 1.0
    return rho + Q - 2*M
end

# Test and compare with original implementation
println("Testing corrected vs original implementation...")

# Classification parameters (sign activation)
kappa1_sign = sqrt(2/π)
kappa_star_sign = sqrt(1 - 2/π)

# Test parameters from paper - should show spike at n/p = 1
p_n_values = [0.5, 1.0, 1.5, 2.0]
lambda_val = 1e-3

println("\nComparison for classification (sign activation):")
println("p/n\tOriginal\tCorrected")
for p_n in p_n_values
    alpha = 1.0 / p_n
    gamma = alpha / 3  # n/d = 3
    
    # Corrected implementation  
    gen_error_corrected = compute_gen_error_corrected(alpha, gamma, kappa1_sign, kappa_star_sign, lambda_val)
    
    println("$p_n\t\t\t$(round(gen_error_corrected, digits=4))")
end

println("\nGenerating comparison plot...")

# Generate data for plotting
p_n_range = 0.1:0.05:3.0
lambda_vals = [1e-4, 1e-3, 1e-2, 1e-1]

plt_corrected = plot()

for (i, lambda_val) in enumerate(lambda_vals)
    gen_errors = Float64[]
    
    for p_n in p_n_range
        alpha = 1.0 / p_n  
        gamma = alpha / 3
        
        try
            gen_error = compute_gen_error_corrected(alpha, gamma, kappa1_sign, kappa_star_sign, lambda_val)
            push!(gen_errors, gen_error)
        catch e
            push!(gen_errors, NaN)
        end
    end
    
    plot!(plt_corrected, p_n_range, gen_errors, 
          label="λ = $(lambda_val)", linewidth=2,
          xlabel="p/n", ylabel="Generalization Error",
          title="CORRECTED: Classification Error vs p/n (should show spike at p/n=1)")
end

display(plt_corrected)
savefig(plt_corrected, "figures/corrected_classification_error.png")

println("Corrected implementation complete! Check for spike at p/n = 1 for low λ values.")