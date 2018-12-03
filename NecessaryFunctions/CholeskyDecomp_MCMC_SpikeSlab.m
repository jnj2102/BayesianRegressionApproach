function [inverse_sigma_squared_vector, indicator_vector_cell,Omega_matrix] = CholeskyDecomp_MCMC_SpikeSlab(Z, p,n, initial_indicator_vector_cell,...
    inverse_sigma_squared_initial, final_rho_cell, sigma_beta_squared,...
    capital_A, capital_B)
%Function to find the sparse precision matrix using a Cholesky
%Decomposition and spike and slab prior
%Author: Jami Jackson Mulgrave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%allocate space

L_offdiag_cell = cell([p-1,1]); 
indicator_vector_cell = cell([p-1,1]); 
inverse_sigma_squared_vector = ones([p,1]);
zeros_accumulation = cell([p-2,1]);

%%%%Sample all offdiagonal columns of the lower triangular matrix

for m = 1:(p-1)  
    %Partition the matrices

Z_exclude_m = Z(:, (m+1):end);
Z_m_vector = Z(:,m);

inverse_sigma_squared = inverse_sigma_squared_initial(m); %inverse sigma squared
rho = final_rho_cell{m};
initial_indicator_vector = initial_indicator_vector_cell{m};

Gamma_diagMatrix = diag(initial_indicator_vector);


%Sample the betas given the sigma_squared and the w_vector

%With the fast sampler algorithm, don't need directly take this inverse.

[numelem, ~] = size(initial_indicator_vector);

IdentityMatrix = eye(numelem);

D_matrix = sigma_beta_squared*IdentityMatrix;

capital_phi = Z_exclude_m*Gamma_diagMatrix*sqrt(inverse_sigma_squared);

alpha = Z_m_vector*sqrt(inverse_sigma_squared);

I_n = eye(n);

%sample u and v
    u = (mvnrnd(zeros([numelem,1]),D_matrix))';

    v = capital_phi*u+ (mvnrnd(zeros([n,1]),I_n))';
    
    %renamed this to l_vector
    l_vector = (capital_phi*D_matrix*capital_phi' + I_n)\(alpha- v);
    
%now sample the beta

    beta_offdiag = u + D_matrix*capital_phi'*l_vector;

%beta_offdiag_cell{m} = beta_offdiag;

% sample the indicators
lambda = log(rho./(1-rho));

indicator_vector = zeros([numelem,1]);

%sample up to (p-1) and sample p separately
for j = 1:numelem
 Z_designMat_removedj = Z_exclude_m;
    Z_designMat_removedj(:,j) = [];
     
    %If these are empty, make these zero (because of the nature of the regression)
    if isempty(Z_designMat_removedj) == 1
        
        Z_designMat_removedj = 0;
    end
     
    indicator_vector_removedj = initial_indicator_vector;
    indicator_vector_removedj(j) = [];
    
        
    if isempty(indicator_vector_removedj) == 1
        
        indicator_vector_removedj = 0;
    end
    
    
    beta_offdiag_removedj = beta_offdiag;
    beta_offdiag_removedj(j) = [];
    
    if isempty(beta_offdiag_removedj) == 1
        
        beta_offdiag_removedj = 0;
    end
    
    
   eta = lambda(j) - 1/2*inverse_sigma_squared*(norm(Z_exclude_m(:,j),2)^2)*...
       beta_offdiag(j)^2   + inverse_sigma_squared*beta_offdiag(j)*...
       Z_exclude_m(:,j)'*(Z_m_vector - Z_designMat_removedj*diag(indicator_vector_removedj)*...
       beta_offdiag_removedj);
   
     Bernoulli_prob = exp(eta)./(1+exp(eta));
     
     %change infinity to one for the indicator.  This can happen sometimes
     %due to the inverse sigma squared value being large.
     
     if exp(eta) == inf
         Bernoulli_prob = 1;
     end
     
      %change really small values to exactly zero
   if Bernoulli_prob < eps
       Bernoulli_prob = 0;
   end
   
          
indicator_vector(j) = binornd(1, Bernoulli_prob);

   
end


indicator_vector_cell{m} = indicator_vector;

%Find the sparse beta
indicatorbeta = beta_offdiag .* indicator_vector;

%finally sample sigma

  %MATLAB uses the scale parameter version of gamma distribution.  To get
  %the inverse, take the reciprocal of the rate parameter,
  %because the sigma^2 ~ InvGamma and 1/(sigma^2) ~ Gamma relationship is
  %dependent on the gamma distribution having a rate parameter.
  
inverse_sigma_squared = gamrnd(capital_A + n/2, 1/(capital_B + 1/2*norm(Z_m_vector -...
    Z_exclude_m*diag(indicator_vector)*beta_offdiag,2)^2));


inverse_sigma_squared_vector(m) = inverse_sigma_squared;

L_diag_exclude_p = sqrt(inverse_sigma_squared);

%Find the L using the most updated beta and inverse_sigma2
L_offdiag = -indicatorbeta.*L_diag_exclude_p;  %Don't need the pth element for the offdiagonals

L_offdiag_cell{m} = L_offdiag;

end %end of d = 1:(p-1)



%Now sample the pth diagonal element

Z_p_vector = Z(:,p);
%there is no dimension p here because we are not sampling beta here and there is no design matrix.  It's
%just Z_p = e_p, where e_p is the error term.
last_inverse_sigma_squared = gamrnd((capital_A + n/2), 1/(capital_B +...
    1/2*norm(Z_p_vector,2)^2));
 %drop the betas here because they are zero

inverse_sigma_squared_vector(p) = last_inverse_sigma_squared;

%Now solve for the diagonal element of lower triangular
%matrix
L_diag_vector = sqrt(inverse_sigma_squared_vector);


%now put together the lower triangular matrix

%make cells of progressively longer zeros
for m = 1:(p-2) 
    
    zeros_accumulation{m} = zeros([m,1]); 
end

%first element is the diagonal element and the offdiagonals.  No zeros.
lower_triangle_first = [L_diag_vector(1);L_offdiag_cell{1}];

lower_triangle_rest = [];

for d = 2:p-1
    
 lower_triangle_rest = horzcat(lower_triangle_rest,[zeros_accumulation{d-1};L_diag_vector(d);L_offdiag_cell{d}]);
    
end

%last element is just the diagonal element and all zeros on top

last_column_zeros = zeros([p-1,1]);
lower_triangle_last_column = [last_column_zeros; L_diag_vector(p)];

lower_triangular_matrix = horzcat(lower_triangle_first,lower_triangle_rest,lower_triangle_last_column);

Omega_matrix = lower_triangular_matrix*lower_triangular_matrix';
    
    

     

    
end















