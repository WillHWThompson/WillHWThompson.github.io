using Plots

# Define the Stieltjes transform g_mu(z) and its derivative
function g_mu(z, γ)
    #sqrt_term = sqrt((z - 1 - γ)^2 - 4 * γ)
    sqrt_term = sqrt((z-1-γ)^2 - 4*γ)
    return (1 - z - γ - sqrt_term) / (2 * z * γ)
end

function g_mu_derivative(z, gamma)
        # Compute the square root term
        sqrt_term = sqrt(-4 * gamma + (1 + gamma - z)^2)
        # Numerator
        numerator = (1 - 2 * gamma + gamma^2 - sqrt_term + gamma * sqrt_term - z - gamma * z)

        denominator = (2 * gamma * sqrt_term * z^2)

        return numerator / denominator
end


function classification_hat_overlaps(alpha,gamma,kappa_1,kappa_star,V_inf,M_inf,Q_inf)
        # Update hatted overlaps
        V_s_hat_inf = (alpha * kappa_1^2) / (gamma * (1 + V_inf))
        q_s_hat_inf = (alpha * kappa_1^2 * (1 + Q_inf - (2 * sqrt(2) * M_inf) / sqrt(π))) / (gamma * (1 + V_inf)^2)
        m_s_hat_inf = (alpha/gamma * sqrt(2 / π) * kappa_1) / (1 + V_inf)
        V_w_hat_inf = (alpha * kappa_star^2) / (1 + V_inf)
        q_w_hat_inf = (alpha * kappa_star^2 * (1 + Q_inf - (2 * M_inf) / sqrt(π))) / (1 + V_inf)^2
        return V_s_hat_inf,q_s_hat_inf,m_s_hat_inf,V_w_hat_inf,q_w_hat_inf
end

function regression_hat_overlaps(alpha,gamma,kappa_1,kappa_star,V_inf,M_inf,Q_inf,Delta)
        # Update hatted overlaps
        V_s_hat_inf = (alpha * kappa_1^2) / (gamma * (1 + V_inf))
        q_s_hat_inf = (alpha * kappa_1^2 * (1 + Delta + Q_inf - 2* M_inf)) / (gamma * (1 + V_inf)^2)
        m_s_hat_inf = (alpha/gamma * kappa_1) / ((1 + V_inf))
        V_w_hat_inf = (alpha * kappa_star^2) / (1 + V_inf)
        q_w_hat_inf = (alpha * kappa_star^2 * (1+ Delta + Q_inf - (2 * M_inf))) / (1 + V_inf)^2
        return V_s_hat_inf,q_s_hat_inf,m_s_hat_inf,V_w_hat_inf,q_w_hat_inf
end


function iterated_equation_classification(alpha, gamma, kappa_1, kappa_star,lambda;Delta = 1,kind = :classification)
    tolerance = 1e-8
    max_iterations = 1000
    # Initialize variables
    V_inf = 0.1
    Q_inf = 0.1
    M_inf = 0.1

    V_s_hat_inf,q_s_hat_inf,m_s_hat_inf,V_w_hat_inf,q_w_hat_inf = 0,0,0,0,0
    for iteration in 1:max_iterations
        # Update hatted overlaps
        # V_s_hat_inf = (alpha * kappa_1^2) / (gamma * (1 + V_inf))
        # q_s_hat_inf = (alpha * kappa_1^2 * (1 + Q_inf - (2 * sqrt(2) * M_inf) / sqrt(π))) / (gamma * (1 + V_inf)^2)
        # m_s_hat_inf = (alpha/gamma * sqrt(2 / π) * kappa_1) / ((1 + V_inf))
        # V_w_hat_inf = (alpha * kappa_star^2) / (1 + V_inf)
        # q_w_hat_inf = (alpha * kappa_star^2 * (1 + Q_inf - (2 * M_inf) / sqrt(π))) / (1 + V_inf)^2

        if kind == :classification
            V_s_hat_inf,q_s_hat_inf,m_s_hat_inf,V_w_hat_inf,q_w_hat_inf = classification_hat_overlaps(alpha,gamma,kappa_1,kappa_star,V_inf,M_inf,Q_inf)
        elseif kind == :regression
            V_s_hat_inf,q_s_hat_inf,m_s_hat_inf,V_w_hat_inf,q_w_hat_inf = regression_hat_overlaps(alpha,gamma,kappa_1,kappa_star,V_inf,M_inf,Q_inf,Delta)
        else 
            println("Error: $kind is not a valid")
        end


        # Compute z for the Stieltjes transform
        z = (lambda + V_w_hat_inf) / V_s_hat_inf 
        g_z = g_mu(-z, gamma)
        g_prime_z = g_mu_derivative(-z, gamma)

        # Compute non-hatted variables (D.23)
        V_s_inf = (1 / V_s_hat_inf) * (1 - z * g_z)#good
        q_s_inf = ((m_s_hat_inf^2 + q_s_hat_inf) / (V_s_hat_inf^2)) * (1 - 2 * z * g_z + z^2 * g_prime_z) -
                  (q_w_hat_inf / ((lambda + V_w_hat_inf) * V_s_hat_inf)) * (-z * g_z + z^2 * g_prime_z)
        m_s_inf = (m_s_hat_inf / V_s_hat_inf) * (1 - z * g_z)

        V_w_inf = gamma / (lambda + V_w_hat_inf) * (1/gamma - 1 + z * g_z)
        q_w_inf = gamma * (q_w_hat_inf/(lambda + V_w_hat_inf)^2) * (1/gamma  - 1 + z^2 * g_prime_z) - gamma*(m_s_hat_inf^2 + q_s_hat_inf)/((lambda + V_w_hat_inf)*V_s_hat_inf)*(-z*g_z + z^2*g_prime_z)
        
        if q_w_inf == Inf
            println("True")
        end

        V_new = kappa_1^2 * V_s_inf + kappa_star^2 * V_w_inf
        Q_new = kappa_1^2 * q_s_inf + kappa_star^2 * q_w_inf
        M_new = kappa_1 * m_s_inf
         


        #println("iteration: $iteration, q_s_inf: $q_s_inf, q_w_inf: $q_w_inf")
        # Check for convergence
        if abs(V_new - V_inf) < tolerance && abs(Q_new - Q_inf) < tolerance && abs(M_new - M_inf) < tolerance
            break
        end

        # Update for the next iteration
        V_inf, Q_inf, M_inf = V_new, Q_new, M_new
    end

    return V_inf, Q_inf, M_inf
end

# Function to compute the generalization error
function compute_generalization_error(alpha,gamma, kappa_1, kappa_star,lambda;kind = :classification,Delta = 0)

    V_inf, Q_inf, M_inf = iterated_equation_classification(alpha, gamma, kappa_1, kappa_star,lambda,kind = kind,Delta = Delta)
    # Compute generalization error
    rho = 1.0  # Assumes ||θ⁰||² / d = 1
    epsilon_g = rho + Q_inf - 2 * M_inf
    return epsilon_g
end




""" Classification W/ ERF
"""
# kappa_1 = sqrt(2 / π)           # Scaling factor for the sign activation
# kappa_star = sqrt(1 - (2 / π))  # Complementary scaling factor


kappa_1_erf = 2/sqrt(3*pi)           # Scaling factor for the sign activation
kappa_star_erf = 0.2003#sqrt(1 - (2 / π))  # Complementary scaling factor

#Parameters
p_n = 0.1:0.01:2  # Range of p/n values
lambda_list = [5e-3,1e-3,1e-2,1,2]
#lambda_list = [5e-3,1e-2,1e-1,1]
results_classification=zeros((size(p_n)[1],size(lambda_list)[1]))
for (i, p_n_val) in enumerate(p_n)
    for (j, lambda_val) in enumerate(lambda_list)
        alpha = 1 / p_n_val
        gamma = alpha/3
        epsilon_list = compute_generalization_error(alpha,gamma ,kappa_1_erf, kappa_star_erf, lambda_val,kind = :classification)
        results_classification[i,j] = epsilon_list
    end
end
# Plot results with labels for lambda values
plt = plot()
colors = cgrad(:bluesreds, length(lambda_list))
for (j, color) in enumerate(colors)
    @show j
    plot!(plt, p_n, results_classification[:, j], label="λ=$(lambda_list[j])", color=color,lw = 2)
end
xlabel!(plt, "p/n")
ylabel!(plt, "Generalization Error")
title!(plt, "Regression Generalization Error vs p/n with Erf(x)")
display(plt)
savefig("figures/classificiation_erf_gen_err.png")



""" Classification W/ sign
"""
kappa_1_sign = sqrt(2 / π)           # Scaling factor for the sign activation
kappa_star_sign = sqrt(1 - (2 / π))  # Complementary scaling factor

#Parameters
p_n = 0.1:0.01:2  # Range of p/n values
lambda_list = [5e-3,1e-3,1e-2,1,2]
#lambda_list = [5e-3,1e-2,1e-1,1]
results_classification=zeros((size(p_n)[1],size(lambda_list)[1]))
for (i, p_n_val) in enumerate(p_n)
    for (j, lambda_val) in enumerate(lambda_list)
        alpha = 1 / p_n_val
        gamma = alpha/3
        epsilon_list = compute_generalization_error(alpha,gamma ,kappa_1_sign, kappa_star_sign, lambda_val,kind = :classification)
        results_classification[i,j] = epsilon_list
    end
end
# Plot results with labels for lambda values
plt = plot()
colors = cgrad(:bluesreds, length(lambda_list))
for (j, color) in enumerate(colors)
    @show j
    plot!(plt, p_n, results_classification[:, j], label="λ=$(lambda_list[j])", color=color,lw = 2)
end
xlabel!(plt, "p/n")
ylabel!(plt, "Generalization Error")
title!(plt, "Regression Generalization Error vs p/n with sgn(x)")
display(plt)
savefig("figures/classificiation_sign_gen_err.png")





"""
Regression with Linear Labels
"""
kappa_0_reg = 0
kappa_1_reg = 1
kappa_star_reg = 0

lambda_list = [5e-4,1e-3,1e-2,1,2]
results_regression=zeros((size(p_n)[1],size(lambda_list)[1]))
for (i, p_n_val) in enumerate(p_n)
    for (j, lambda_val) in enumerate(lambda_list)
        alpha = 1 / p_n_val
        gamma = 0.8#alpha/
        epsilon_list = compute_generalization_error(alpha,gamma ,kappa_1_reg, kappa_star_reg, lambda_val,kind = :regression,Delta=0.1)
        results_regression[i,j] = epsilon_list
    end
end

plt = plot()
colors = cgrad(:bluesreds, length(lambda_list))
for (j, color) in enumerate(colors)
    plot!(plt, p_n, results_regression[:, j], label="λ=$(lambda_list[j])", color=color,lw = 2)
end
xlabel!(plt, "p/n")
ylabel!(plt, "Generalization Error")
title!(plt, "Generalization Error vs p/n")
display(plt)

savefig("figures/regression_gen_err.png")


"""
Regression with Linear Labels
"""
kappa_0_reg = 0
kappa_1_reg = 1
kappa_star_reg = 0

Delta_list = [1e-2,1e-1,1]
results_regression=zeros((size(p_n)[1],size(Delta_list)[1]))
for (i, p_n_val) in enumerate(p_n)
    for (j, Delta_val) in enumerate(Delta_list)
        alpha = 1 / p_n_val
        gamma = 0.8#alpha/
        epsilon_list = compute_generalization_error(alpha,gamma ,kappa_1_reg, kappa_star_reg, 1e-1,kind = :regression,Delta=Delta_val)
        results_regression[i,j] = epsilon_list
    end
end

plt = plot()
colors = cgrad(:bluesreds, length(Delta_list))
for (j, color) in enumerate(colors)
    plot!(plt, p_n, results_regression[:, j], label="Δ=$(Delta_list[j])", color=color,lw = 2)
end
xlabel!(plt, "p/n")
ylabel!(plt, "Generalization Error")
title!(plt, "Generalization Error vs p/n")
display(plt)

savefig("figures/regression_gen_err_noise_sweep.png")








p_n = 0.1:0.01:2   # Range of p/n values
n_d = 0.5:0.01:2   # Range of n/d values
lambda_val = 1e-1  # Lambda value
# Initialize the results matrix
# Note: Swap the dimensions to match the new axes
results_nd_pn_heatmap = zeros(length(n_d), length(p_n))

lambda_val = 1e-2
# Swap the loops: n_d is now the outer loop
for (i, n_d_val) in enumerate(n_d)
    for (j, p_n_val) in enumerate(p_n)
        # Compute alpha and gamma
        alpha = 1.0 / p_n_val
        gamma =  alpha/ n_d_val
        # Compute the generalization error
        epsilon = compute_generalization_error(alpha, gamma, kappa_1_sign, kappa_star_erf, lambda_val)
        # Store the result
        results_nd_pn_heatmap[i, j] = epsilon
    end
end

# Create the surface plot with flipped axes
surface_plot = surface(
    n_d,   # n/d on the x-axis
    p_n,   # p/n on the y-axis
    results_nd_pn_heatmap',  # Use the results matrix as is
    xlabel = "n/d",
    ylabel = "p/n",
    zlabel = "Generalization Error",
    title = "Generalization Error Classification with Erf"
)

# Adjust the camera angle if desired
plot!(surface_plot, camera = (80, 30))
# Display the plot
display(surface_plot)

savefig("figures/surface_plot_nd_pn_erf_classification.png")





# #Compute generalization error for different p/n ratiosk


# #
# epsilon_list = compute_generalization_error(p_n, kappa_1, kappa_star,lambda)
# ymin = 0
# ymax=1

# # Plot generalization error
# plt = plot(
#     p_n, 
#     epsilon_list, 
#     xlabel="p/n", 
#     ylabel="Generalization Error", 
#     title="Generalization Error vs p/n", 
#     ylim=(ymin, ymax)
# )
# display(plt)
# println("done")

# alpha_i = 1.0
# gamma_i = 0.1

# kappa_1 = sqrt(2 / π)
# kappa_star = sqrt(1 - (2 / π))
# lambda = 1.4
# rho = 1

# epsilon_g_list = []
# alpha_vals = 0.1:0.01:2
# for alpha_i in alpha_vals
#     V_inf, Q_inf, M_inf = iterated_equation_classification(alpha_i, gamma_i, kappa_1, kappa_star, lambda)
#     epsilon_g = rho + Q_inf - 2 * M_inf
#     push!(epsilon_g_list,epsilon_g)
# end

# plot(alpha_vals,epsilon_g_list,title = "test")


#println("V_inf = $V_inf, Q_inf = $Q_inf, M_inf = $M_inf