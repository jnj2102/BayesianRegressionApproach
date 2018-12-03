function [lower_bound_temp] = LowerboundMarginalLikelihoodCholeskyDecomp(Z_designMat,...
    Z_outcome, zeta_squared, capital_A, capital_B, tau, w_vector, n, rho)
    
% Calculate the lower bound of the marginal likelihood for the variational
% Bayes algorithm to do spike and slab with the Cholesky Decomposition
% Author: Jami Jackson Mulgrave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
W_diagMatrix = diag(w_vector);

[dimension, ~] = size(W_diagMatrix);

XtX = (Z_designMat'*Z_designMat);

Identity_matrix = eye(dimension);

D_diagMatrix = tau*XtX .* W_diagMatrix .* (Identity_matrix - W_diagMatrix) + ...
    1/zeta_squared*Identity_matrix;

Sigma_beta = inv(tau*W_diagMatrix*XtX*W_diagMatrix + D_diagMatrix);

mean_beta = tau*Sigma_beta*W_diagMatrix*Z_designMat'*Z_outcome;

Omega_matrix = w_vector*w_vector' + W_diagMatrix*(Identity_matrix - W_diagMatrix);

s_inverseGamma = capital_B + 1/2*(norm(Z_outcome,2)^2 - 2*Z_outcome'*...
    Z_designMat*W_diagMatrix*mean_beta + trace((XtX .* Omega_matrix)*...
    (mean_beta*mean_beta' + Sigma_beta)));


%Check if each element of w vector is on the boundary

element_Bernoulli = zeros([dimension, 1]);

for m = 1:dimension
    if w_vector(m) == 0
        
    element_Bernoulli(m) = (1-w_vector(m))*log((1-rho(m))/...
    (1-w_vector(m)));

    elseif w_vector(m) == 1
            
      element_Bernoulli(m) =  w_vector(m)*log(rho(m)/w_vector(m));
      
    else
        
      element_Bernoulli(m) = w_vector(m)*log(rho(m)/w_vector(m)) +...
          (1-w_vector(m))*log((1-rho(m))/(1-w_vector(m)));
    end
    
end
    
sum_Bernoulli = sum(element_Bernoulli);

logdet_Sigma_beta = 2 * sum(log(diag(chol(Sigma_beta))));

lower_bound_temp = dimension/2 -n/2*log(2*pi) +...
    -(dimension/2)*log(zeta_squared) +...
 + capital_A*log(capital_B) +...
    -log(gamma(capital_A)) + log(gamma(capital_A + n/2)) +...
    -(capital_A + n/2)*log(s_inverseGamma) + 1/2*logdet_Sigma_beta + ...
    -1/(2*zeta_squared)*trace(mean_beta*mean_beta' + Sigma_beta) +...
    sum_Bernoulli;
    

    

    
    
end


