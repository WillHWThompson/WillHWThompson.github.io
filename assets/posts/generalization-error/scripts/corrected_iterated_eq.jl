using Plots
using QuadGK
using Distributions

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

# Define the functions η(y,ω) and Z(y,ω) for different loss functions
# Classification with square loss and sign labels
function eta_classification_square(y, omega, V_val)
    # For classification: y ∈ {-1, +1}, loss is (y - x)²/2
    # η(y,ω) = argmin_x [(x-ω)²/(2V) + (y-x)²/2]
    return (omega/V_val + y) / (1/V_val + 1)
end

function Z_classification_square(y, omega, V_val, rho_val, M_val, Q_val)
    # For sign labels: y = sign(cμ·θ⁰/√d) 
    # Z(y,ω) integrates over the label distribution
    V_prime = rho_val - M_val^2/Q_val
    if V_prime <= 0
        return 1e-10  # Avoid numerical issues
    end
    
    omega_prime = M_val/sqrt(Q_val)
    
    # For sign labels, we integrate the Gaussian distribution
    # P_y^0(y|ν) for ν ~ N(0, V_prime)
    result, _ = quadgk(nu -> begin
        prob_nu = exp(-nu^2 / (2*V_prime)) / sqrt(2π * V_prime)
        if y * nu > 0  # Sign matches
            prob_nu
        else
            0.0
        end
    end, -Inf, Inf, rtol=1e-6)
    
    return max(result, 1e-10)  # Avoid numerical issues
end

# Regression with square loss and linear labels  
function eta_regression_square(y, omega, V_val)
    # For regression: η(y,ω) = argmin_x [(x-ω)²/(2V) + (y-x)²/2]
    return (omega/V_val + y) / (1/V_val + 1)
end

function Z_regression_square(y, omega, V_val, rho_val, M_val, Q_val, Delta=0.0)
    # For linear labels with noise: y = cμ·θ⁰/√d + √Δ ξ
    V_prime = rho_val - M_val^2/Q_val + Delta
    if V_prime <= 0
        return 1e-10
    end
    
    omega_prime = M_val/sqrt(Q_val)
    
    # For linear labels: P_y^0(y|ν) = N(y; ν, Δ)
    prob = exp(-(y - omega_prime)^2 / (2*V_prime)) / sqrt(2π * V_prime)
    return max(prob, 1e-10)
end

# Fixed point iteration implementing equations D.22-D.25
function fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda; 
                              loss_type=:classification_square, Delta=0.0,
                              max_iter=1000, tol=1e-8)
    
    # Initialize
    V_s_hat, q_s_hat, m_s_hat = 0.1, 0.1, 0.1
    V_w_hat, q_w_hat = 0.1, 0.1
    
    rho = 1.0  # ||θ⁰||²/d
    
    for iter in 1:max_iter
        # Compute V, Q, M from hat variables
        V_val = kappa1^2 * (1/V_s_hat * (1 - (lambda + V_w_hat)/V_s_hat * g_mu(-(lambda + V_w_hat)/V_s_hat, gamma))) + 
                kappa_star^2 * (gamma/(lambda + V_w_hat) * (1/gamma - 1 + (lambda + V_w_hat)/V_s_hat * g_mu(-(lambda + V_w_hat)/V_s_hat, gamma)))
        
        z = (lambda + V_w_hat) / V_s_hat
        g_z = g_mu(-z, gamma)
        g_prime_z = g_mu_derivative(-z, gamma)
        
        V_s = (1/V_s_hat) * (1 - z * g_z)
        q_s = ((m_s_hat^2 + q_s_hat)/(V_s_hat^2)) * (1 - 2*z*g_z + z^2*g_prime_z) - 
              (q_w_hat/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
        m_s = (m_s_hat/V_s_hat) * (1 - z * g_z)
        
        V_w = gamma/(lambda + V_w_hat) * (1/gamma - 1 + z*g_z)
        q_w = gamma * (q_w_hat/(lambda + V_w_hat)^2) * (1/gamma - 1 + z^2*g_prime_z) + 
              ((m_s_hat^2 + q_s_hat)/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
        
        Q_val = kappa1^2 * q_s + kappa_star^2 * q_w
        M_val = kappa1 * m_s
        
        # Numerical integration for the hat variables (D.22)
        # Sample ξ ~ N(0,1) for Monte Carlo integration
        n_samples = 1000
        samples_xi = randn(n_samples)
        
        # For y sampling, we need to be careful about the label distribution
        if loss_type == :classification_square
            y_values = [-1, 1]  # Binary classification
        else
            y_values = randn(100)  # Sample from continuous distribution for regression
        end
        
        # Compute integrals via Monte Carlo
        V_s_hat_new = 0.0
        q_s_hat_new = 0.0
        m_s_hat_new = 0.0
        V_w_hat_new = 0.0
        q_w_hat_new = 0.0
        
        for xi in samples_xi
            omega0 = M_val/sqrt(Q_val) * xi
            omega1 = sqrt(Q_val) * xi
            
            for y in y_values
                if loss_type == :classification_square
                    Z_val = Z_classification_square(y, omega0, V_val, rho, M_val, Q_val)
                    eta_val = eta_classification_square(y, omega1, V_val)
                    weight = 0.5  # Equal probability for y = ±1
                else
                    Z_val = Z_regression_square(y, omega0, V_val, rho, M_val, Q_val, Delta)
                    eta_val = eta_regression_square(y, omega1, V_val)
                    weight = exp(-y^2/2) / sqrt(2π)  # Weight for sampling y
                end
                
                # Derivatives for D.22
                deta_domega = 1 / (1/V_val + 1)  # For square loss
                
                V_s_hat_new += (alpha * kappa1^2 / gamma) * Z_val * deta_domega * weight
                q_s_hat_new += (alpha * kappa1^2 / gamma) * Z_val * ((eta_val - omega1)/V_val)^2 * weight
                m_s_hat_new += (alpha * kappa1 / gamma) * Z_val * ((eta_val - omega1)/V_val) * (omega0/sqrt(rho)) * weight
                V_w_hat_new += alpha * kappa_star^2 * Z_val * deta_domega * weight
                q_w_hat_new += alpha * kappa_star^2 * Z_val * ((eta_val - omega1)/V_val)^2 * weight
            end
        end
        
        V_s_hat_new /= n_samples * length(y_values)
        q_s_hat_new /= n_samples * length(y_values)
        m_s_hat_new /= n_samples * length(y_values)
        V_w_hat_new /= n_samples * length(y_values)
        q_w_hat_new /= n_samples * length(y_values)
        
        # Check convergence
        if (abs(V_s_hat_new - V_s_hat) < tol && abs(q_s_hat_new - q_s_hat) < tol && 
            abs(m_s_hat_new - m_s_hat) < tol && abs(V_w_hat_new - V_w_hat) < tol && 
            abs(q_w_hat_new - q_w_hat) < tol)
            break
        end
        
        # Update with damping for stability
        damping = 0.1
        V_s_hat = (1-damping) * V_s_hat + damping * V_s_hat_new
        q_s_hat = (1-damping) * q_s_hat + damping * q_s_hat_new  
        m_s_hat = (1-damping) * m_s_hat + damping * m_s_hat_new
        V_w_hat = (1-damping) * V_w_hat + damping * V_w_hat_new
        q_w_hat = (1-damping) * q_w_hat + damping * q_w_hat_new
    end
    
    # Compute final values
    z = (lambda + V_w_hat) / V_s_hat
    g_z = g_mu(-z, gamma)
    g_prime_z = g_mu_derivative(-z, gamma)
    
    V_s = (1/V_s_hat) * (1 - z * g_z)
    q_s = ((m_s_hat^2 + q_s_hat)/(V_s_hat^2)) * (1 - 2*z*g_z + z^2*g_prime_z) - 
          (q_w_hat/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
    m_s = (m_s_hat/V_s_hat) * (1 - z * g_z)
    
    V_w = gamma/(lambda + V_w_hat) * (1/gamma - 1 + z*g_z)
    q_w = gamma * (q_w_hat/(lambda + V_w_hat)^2) * (1/gamma - 1 + z^2*g_prime_z) + 
          ((m_s_hat^2 + q_s_hat)/((lambda + V_w_hat)*V_s_hat)) * (-z*g_z + z^2*g_prime_z)
    
    Q_val = kappa1^2 * q_s + kappa_star^2 * q_w
    M_val = kappa1 * m_s
    V_final = kappa1^2 * V_s + kappa_star^2 * V_w
    
    return V_final, Q_val, M_val
end

# Compute generalization error
function compute_gen_error_corrected(alpha, gamma, kappa1, kappa_star, lambda; 
                                   loss_type=:classification_square, Delta=0.0)
    V, Q, M = fixed_point_iteration(alpha, gamma, kappa1, kappa_star, lambda; 
                                   loss_type=loss_type, Delta=Delta)
    rho = 1.0
    return rho + Q - 2*M
end

# Test the corrected implementation
println("Testing corrected implementation...")

# Classification parameters (sign activation)
kappa1_sign = sqrt(2/π)
kappa_star_sign = sqrt(1 - 2/π)

# Test single point
alpha_test = 1.0  # n/p = 1
gamma_test = alpha_test / 3  # n/d = 3  
lambda_test = 1e-3

println("Testing single point: α=$alpha_test, γ=$gamma_test, λ=$lambda_test")

try
    gen_error = compute_gen_error_corrected(alpha_test, gamma_test, kappa1_sign, kappa_star_sign, lambda_test)
    println("Generalization error: $gen_error")
catch e
    println("Error: $e")
end

println("Corrected implementation created!")