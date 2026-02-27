using Random, Distributions
using LinearAlgebra
using Flux
using Statistics  # For mean and std

# Add the src directory to the load path
push!(LOAD_PATH, "../src")


using Flux
using Flux.Optimise: update!, Adam  # Add Adam to the imports
using Flux: gradient, params
using Infiltrator


"""
Code for Data Generation
"""
# Include and use the module
function haar_orthogonal_matrix(n, p)
    # sample uniformly from the space of orthogonal matrices (the Haar measure)
    A = randn(n, p)
    Q, R = qr(A)
    Q = Q ./ abs.(diag(R))
    return Q
end

function sample_F(n, p, d, γ; method = :uniform)
    if method == :rand_orthogonal_projections # sample matrix via random orthogonal projection
        U = haar_orthogonal_matrix(d, d)
        V = haar_orthogonal_matrix(p, p)
        # generate D
        d_k = max(sqrt(γ), 1)
        D = d_k .* Matrix{Float64}(I, d, p)
        F = U * D * V
        return F
    elseif method == :uniform
        return randn(d, p)
    end
end


function generate_dataset(d, α, γ, f_0, my_σ,θ₀,F)
    p = floor(Int, d / γ)
    n = floor(Int, α * p)

    A = randn(n, p)
    C = randn(n, d)

    y = f_0(C * θ₀') # generate a matrix of n data points
    X = my_σ(C * F)
    return X, y
end


# Define the GLM model
mutable struct GLM
    weights::Vector
    link_function::Function  # The link function g()
    inverse_link_function::Function  # The inverse link function g⁻¹()
end

# Define the prediction function for the GLM
function predict(glm::GLM, X)
    glm.inverse_link_function.(X * glm.weights)
end

# Define the loss function
function loss_function(glm::GLM, X, y, λ)
    # Predictions using the model
    predictions = predict(glm, X)
    
    mse_loss = sum((y .- predictions).^2) / (2 * size(X, 1))  # Mean squared error
    regularization = λ * sum(glm.weights.^2) / 2  # Ridge regularization
    return mse_loss + regularization
end

# Define the logistic loss function for binary classification
function logistic_loss(glm::GLM, X, y, λ)
    z = X * glm.weights
    log_likelihood = sum(log(1 .+ exp.(-y .* z))) / size(X, 1)
    regularization = λ * sum(glm.weights.^2) / 2
    return log_likelihood + regularization
end

# Train the GLM using gradient descent
function train_glm!(glm::GLM, X, y, loss_fn, λ; epochs=100000, lr=0.1, tol=1e-4)
    opt = Adam(lr)  # Use Adam optimizer instead of Descent
    ps = params(glm.weights)  # Collect parameters
    
    for epoch in 1:epochs
        # Compute the gradient of the loss with respect to the weights
        gs = gradient(ps) do
            loss_fn(glm, X, y, λ)
        end
        # Update the weights using the optimizer
        update!(opt, ps, gs)
        
        # Check for convergence
        grad_norm = maximum(abs.(gs[glm.weights]))
        if grad_norm < tol
            println("Converged after $epoch epochs")

            my_loss = loss_fn(glm, X, y, λ)
            println("Loss of $my_loss")
            break
        end

    if epoch == epochs
        println("Did not converge")
    end
    end
end

# Function to calculate R^2
function r_squared(y_true, y_pred)
    ss_total = sum((y_true .- mean(y_true)).^2)
    ss_res = sum((y_true .- y_pred).^2)
    return 1 - (ss_res / ss_total)
end

function generlization_error(glm,d,alpha,gamma,f_0,my_σ;task = :regression)
    #calculate the generalization error
    if task == :regression
        k = 0
    elseif  task == :classification
        k = 1
    end

    X_new,y_new = generate_dataset(d,alpha,gamma,identity,identity)
    gen_error_unnormed = loss_function(glm,X_new,y_new,λ)
    return 1/(n*4^k )* gen_error_unnormed
end

function averaged_generalization_error(d, α, γ, λ, my_σ, f₀, f̂, nseeds; task = :regression)

    p = floor(Int, d / γ)
    n = floor(Int, α * p)
    Eg = 0.0

    # @show p
    # @show n

    for i in 1:nseeds
        @show i

        θ₀ = randn(1, d)
        F = sample_F(n, p, d, γ, method = :uniform) # sample a matrix for the lifting operator F, can be iid Gaussian or a random orthogonal projection


        # Generate datasets X, y
        X, y = generate_dataset(d, α, γ, f₀, my_σ,θ₀,F)

        @show size(X)

        # Initialize and train the GLM to compute ŵ
        glm = GLM(zeros(p), x -> x, x -> x)

        train_glm!(glm, X, y, loss_function, λ; epochs = 10000, lr = 0.01)

        #@infiltrate

        # Generate new datasets X_new, y_new
        X_new, y_new = generate_dataset(d, α, γ, f₀, my_σ,θ₀,F)

        # Compute predictions on new data
        y_pred = predict(glm, X_new)

        r2 = r_squared(y_new, y_pred)
        println("R^2: $r2")

        #@show n
        # Compute and accumulate the error
        if task == :regression
            k = 0
            error = sum((y_new .- y_pred).^2) / (4^k * n)
            @show error
        elseif task == :classification
            k = 1
            error = sum(log(1 .+ exp.(-y_new .* (X_new * glm.weights)))) / (4^k * n)
        end
        Eg += error
    end

    # Return the averaged generalization error

    εg = Eg / nseeds

    #return εg
    return εg
end

# d = 10
# alpha = 10
# gamma = 1
# lambda = 1e-4



# n_p = 0.1:1:5
# n_seeds = 2


#gen_error_list = []
# for n_p_i in n_p
#     @show n_p_i
#     alpha_i = 1/n_p_i
#     gamma_i = alpha ./ 3
#     epsilon_g =  averaged_generalization_error(d,alpha,gamma,lambda,identity,identity,identity,n_seeds)
#     push!(gen_error_list,epsilon_g)
# end

d = 10
alpha_i = 20
gamma_i = 1
lambda = 1e-4
n_seeds = 500

# epsilon_g =  averaged_generalization_error(d,alpha,gamma,lambda,identity,identity,identity,n_seeds)
# @show epsilon_g
# #push!(gen_err_list,epsilon_g)

gen_err_list = []
alpha_list = 1:1:50
for alpha_i in alpha_list
    epsilon_g =  averaged_generalization_error(d,alpha,gamma,lambda,identity,identity,identity,n_seeds)
    @show epsilon_g
    push!(gen_err_list,epsilon_g)
end


# # # using Plots
plot(alpha_list,gen_err_list)
#push!(gen_error_list,epsilon_g)




γ = gamma 
α = alpha_i
λ=1e-4
my_σ = identity

p = floor(Int, d / γ)
n = floor(Int, α * p)

f₀ = identity
#my_σ = identity




θ₀ = randn(1, d)
F = sample_F(n, p, d, γ, method = :uniform) # sample a matrix for the lifting operator F, can be iid Gaussian or a random orthogonal projection


# Generate datasets X, y
X, y = generate_dataset(d, α, γ, f₀, my_σ,θ₀,F)

@show size(X)

# Initialize and train the GLM to compute ŵ
glm = GLM(zeros(p), x -> x, x -> x)

train_glm!(glm, X, y, loss_function, λ; epochs = 1000000, lr = 0.01)

#@infiltrate

# Generate new datasets X_new, y_new
X_new, y_new = generate_dataset(d, α, γ, f₀, my_σ,θ₀,F)

# Compute predictions on new data
y_pred = predict(glm, X_new)

k = 0
my_error = sum((y_new .- y_pred).^2) / (4^k * n)


r2 = r_squared(y_new, y_pred)
println("R^2: $r2")
println("error: $my_error")
