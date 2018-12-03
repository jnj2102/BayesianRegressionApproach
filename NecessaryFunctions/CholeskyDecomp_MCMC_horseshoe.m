function [inverse_lambda_squared_vector,inverse_a_offdiag_vector,b_offdiag_cell,inverse_c_offdiag_cell,...
    inverse_sigma_squared_vector,Omega_matrix] = CholeskyDecomp_MCMC_horseshoe(Z,p,...
    inverse_lambda_squared_initial, inverse_sigma_squared_initial,...
    b_offdiag_initial,inverse_a_offdiag_initial, capital_A,inverse_c_offdiag_initial, a0, b0,n,...
    c_tune)

%Function to run the full MCMC version of Cholesky decomposition with
%horseshoe prior on regression coefficients
%Author: Jami Jackson Mulgrave


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%allocate space

L_offdiag_cell = cell([p-1,1]); 
inverse_lambda_squared_vector = ones([p-1,1]);
b_offdiag_cell = cell([p-1,1]);
inverse_c_offdiag_cell = cell([p-1,1]);
beta_offdiag_cell = cell([p-1,1]);
inverse_sigma_squared_vector = ones([p,1]);
inverse_a_offdiag_vector = ones([p-1,1]);
zeros_accumulation = cell([p-2,1]);


%%%%Sample all offdiagonal columns of the lower triangular matrix

for m = 1:(p-1)  
%Partition the matrices

Z_exclude_m = Z(:, (m+1):end);
Z_m_vector = Z(:,m);

inverse_lambda_squared = inverse_lambda_squared_initial(m);
lambda_squared = 1/inverse_lambda_squared;

b_offdiag = b_offdiag_initial{m};

inverse_c_offdiag = inverse_c_offdiag_initial{m};

inverse_a_offdiag = inverse_a_offdiag_initial(m);

row_offdiag = (m+1:p)';  %Need the loop indexed

inverse_sigma_squared = inverse_sigma_squared_initial(m); %inverse sigma squared

diag_matrix = diag((lambda_squared*b_offdiag*c_tune^2)./(p^2*row_offdiag));

%%(1) Sample the betas given sigma

%With the fast sampler algorithm, don't need directly take this inverse.

D_matrix = 1/inverse_sigma_squared*diag_matrix;
capital_phi = sqrt(inverse_sigma_squared)*Z_exclude_m;
alpha = sqrt(inverse_sigma_squared)*Z_m_vector;
I_n = eye(n);

[numelem, ~] = size(D_matrix);

%sample u and v
    u= (mvnrnd(zeros([numelem,1]),D_matrix))';

    v=capital_phi*u+ (mvnrnd(zeros([n,1]),I_n))';

    w_vector = (capital_phi*D_matrix*capital_phi' + I_n)\(alpha- v);
    
%now sample the beta

    beta_offdiag = u + D_matrix*capital_phi'*w_vector;

beta_offdiag_cell{m} = beta_offdiag;



%given beta, sample lambda

  %MATLAB uses the scale parameter version of gamma distribution.  To get
  %the inverse, take the reciprocal of the rate parameter,
  %because the sigma^2 ~ InvGamma and 1/(sigma^2) ~ Gamma relationship is
  %dependent on the gamma distribution having a rate parameter.

inverse_diag_matrix_Forlambda = diag(inverse_sigma_squared*...
    ((p^2*row_offdiag)./(b_offdiag*c_tune^2)));

inverse_lambda_squared = gamrnd(numelem/2 + 1/2, 1./(0.5*beta_offdiag'*inverse_diag_matrix_Forlambda*...
    beta_offdiag + inverse_a_offdiag));

inverse_lambda_squared_vector(m) = inverse_lambda_squared;

%sample a

inverse_a_offdiag = gamrnd(1, (1/(inverse_lambda_squared + capital_A^2)));

inverse_a_offdiag_vector(m) = inverse_a_offdiag;

%sample b

inverse_b_offdiag = gamrnd(1,...
    1./(0.5*inverse_sigma_squared*inverse_lambda_squared*p^2*1/(c_tune^2)*...
    row_offdiag.*beta_offdiag.^2 +   inverse_c_offdiag));



b_offdiag = 1./inverse_b_offdiag;

b_offdiag_cell{m} = b_offdiag;

%sample c

inverse_c_offdiag = gamrnd(1, 1./(inverse_b_offdiag + 1));

inverse_c_offdiag_cell{m} = inverse_c_offdiag;

%finally sample sigma squared
inverse_diag_matrix_ForSigma = diag(inverse_lambda_squared*p^2*1/(c_tune^2)*...
    row_offdiag.*inverse_b_offdiag);

inverse_sigma_squared = gamrnd((a0 + (n+ numelem)/2),...
    1/(b0 + 0.5*beta_offdiag'*(inverse_diag_matrix_ForSigma*beta_offdiag) +...
     0.5*norm((Z_m_vector - Z_exclude_m*beta_offdiag),2)^2));
 
 
inverse_sigma_squared_vector(m) = inverse_sigma_squared;

L_diag_exclude_p = sqrt(inverse_sigma_squared);


%Find the L using the most updated beta and inverse_sigma2
L_offdiag = -beta_offdiag.*L_diag_exclude_p;  %Don't need the pth element for the offdiagonals

L_offdiag_cell{m} = L_offdiag;

end


%Now sample the pth diagonal element

Z_p_vector = Z(:,p);
%there is no dimension p here because we are not sampling beta here and there is no design matrix.  It's
%just Z_p = e_p, where e_p is the error term.
last_inverse_sigma_squared = gamrnd((n/2+ a0), 1/(b0 + 1/2*norm(Z_p_vector,2)^2));
 %drop the betas here because they are zero

inverse_sigma_squared_vector(p) = last_inverse_sigma_squared;

%Now solve for the final diagonal element of lower triangular
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
