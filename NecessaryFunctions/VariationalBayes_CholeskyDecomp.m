function [final_w_vector_d,s_inverseGamma_d, s_inverseGamma_p, mean_beta_d,...
    Sigma_beta_d] = VariationalBayes_CholeskyDecomp(Z, zeta_squared,...
    capital_A, capital_B, tau, final_w_vector_cell, n, p,final_rho_cell)
% Function to do Variational Bayes spike and slab over Cholesky
% Decomposition
% Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


max_while = 100;

epsilon = 1e-6;

tau_vector = tau*ones([p-1,1]);

w_vector_star_cell = cell([p-1,1]);
Sigma_beta_d = cell([p-1,1]);
mean_beta_d =  cell([p-1,1]);
s_inverseGamma_d = zeros([p-1,1]);
lower_bound_vector_d = zeros([p-1,1]);

t = 1;


while (t <= max_while )
    
    for d = 1: (p-1)
        
         Z_designMat = Z(:, (d+1):end);
        Z_outcome = Z(:,d);
        rho = final_rho_cell{d};
        tau_d = tau_vector(d);
        
        lambda = log(rho./(1-rho));

    w_vector = final_w_vector_cell{d};
    
W_diagMatrix = diag(w_vector);

[dimension, ~] = size(w_vector);

Identity_matrix = eye(dimension);

Omega_matrix = w_vector*w_vector' + W_diagMatrix*(Identity_matrix - W_diagMatrix);

XtX = (Z_designMat'*Z_designMat);

D_diagMatrix = tau_d*XtX .* W_diagMatrix .* (Identity_matrix - W_diagMatrix) + ...
    1/zeta_squared*Identity_matrix;

Sigma_beta = inv(tau_d*W_diagMatrix*XtX*W_diagMatrix + D_diagMatrix);

Sigma_beta_d{d} = Sigma_beta;

XtY = Z_designMat'*Z_outcome;

mean_beta = tau_d*Sigma_beta*W_diagMatrix*XtY;

mean_beta_d{d}= mean_beta;

s_inverseGamma = capital_B + 1/2*(norm(Z_outcome,2)^2 - 2*XtY'*W_diagMatrix*mean_beta + trace((XtX .* Omega_matrix)*...
    (mean_beta*mean_beta' + Sigma_beta)));

s_inverseGamma_d(d) = s_inverseGamma;

tau_d = (capital_A + n/2)/s_inverseGamma;

w_vector_star = w_vector;



for j = 1:dimension
    Z_designMat_removedj = Z_designMat;
    Z_designMat_removedj(:,j) = [];
    
    %If these are empty, make these zero (needed for numerical reasons)
    if isempty(Z_designMat_removedj) == 1
        
        Z_designMat_removedj = 0;
    end
    
    w_vector_star_removedj = w_vector_star;
    w_vector_star_removedj(j) = [];
    
        
    if isempty(w_vector_star_removedj) == 1
        
        w_vector_star_removedj = 0;
    end
    
    
    mean_beta_removedj = mean_beta;
    mean_beta_removedj(j) = [];
    
    if isempty(mean_beta_removedj) == 1
        
        mean_beta_removedj = 0;
    end
    
    Sigma_beta_removedj = Sigma_beta;
    Sigma_beta_removedj(j,:) = [];
    
    if isempty(Sigma_beta_removedj) == 1
        
        Sigma_beta_removedj = 0;
    end
    
   eta = lambda(j) - 1/2*tau_d*((mean_beta(j))^2 + Sigma_beta(j,j))*...
       (norm(Z_designMat(:,j),2)^2) + tau_d*Z_designMat(:,j)'*(Z_outcome*...
       mean_beta(j) - Z_designMat_removedj*diag(w_vector_star_removedj)*...
       (mean_beta_removedj*mean_beta(j) + Sigma_beta_removedj(:,j)));
   
   
       
   w_vector_star_temp = exp(eta)./(1+exp(eta));
   
   %change really small values to exactly zero and really large values to 1
   %for numerical stability.
   if w_vector_star_temp < eps
       w_vector_star_temp = 0;
   elseif exp(eta) == inf
       w_vector_star_temp = 1;
   end       
   
   
   w_vector_star(j) = w_vector_star_temp;
end



%Calculate the lower bound and compare for each regression
lower_bound_vector_d(d) = LowerboundVB(Sigma_beta,mean_beta,s_inverseGamma,...
    zeta_squared, capital_A, capital_B, w_vector,n, rho);

w_vector_star_cell{d} = w_vector_star;

tau_vector(d) = tau_d; %put the new dth tau in the vector to use again.

    end %end of d=1:p-1

%Find the pth regression and sigma_p^2
  Z_outcome_p = Z(:,p);
s_inverseGamma_p = capital_B + 1/2*norm(Z_outcome_p,2)^2 ;


lower_bound_p = - n/2*log(2*pi) + capital_A*log(capital_B) -...
    log(gamma(capital_A)) + log(gamma(capital_A + n/2)) -...
    (capital_A + n/2)*log(s_inverseGamma_p);

%Finally sum the lower bounds for the p regressions

final_lower_bound_sum(t) = sum(lower_bound_vector_d) + lower_bound_p;


if t == 1
   final_w_vector_cell = w_vector_star_cell;
    t = t+1;                
    continue;
    
elseif abs(final_lower_bound_sum(t) - final_lower_bound_sum(t-1)) < epsilon
  final_w_vector_d = final_w_vector_cell;
    break;
else
   final_w_vector_cell = w_vector_star_cell;
    t = t+1;     
    continue;
end

end %end while loop

if t == max_while
    fprintf('The algorithm did not converge');

end
    
   






end














