function [lower_bound_temp] = LowerboundVB(Sigma_beta,mean_beta,s_inverseGamma,...
     zeta_squared, capital_A, capital_B, w_vector, n, rho)
    
% Calculate the lower bound of the marginal likelihood for the variational
% Bayes algorithm to do spike and slab with the Cholesky Decomposition
%This is for the algorithm and not for tuning.
% Author: Jami Jackson Mulgrave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dimension, ~] = size(w_vector);

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


