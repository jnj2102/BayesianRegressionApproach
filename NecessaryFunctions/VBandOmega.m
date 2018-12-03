function[Omega_matrix_VB_samples] = VBandOmega(Z, zeta_squared,...
    capital_A, capital_B, tau, final_w_vector_cell, n, p,final_rho_cell, num_VB_samples)


[final_w_vector_d,s_inverseGamma_d, s_inverseGamma_p, mean_beta_d,...
    Sigma_beta_d] = VariationalBayes_CholeskyDecomp(Z, zeta_squared,...
    capital_A, capital_B, tau, final_w_vector_cell, n, p,final_rho_cell);


Omega_matrix_VB_samples = zeros([p,p,num_VB_samples]);

for samps = 1:num_VB_samples
%Construct the precision matrix from the lower triangular matrix

inverse_sigma_squared_vector = zeros([p,1]);

L_offdiag_cell = cell([p-1,1]); 

zeros_accumulation = cell([p-2,1]);

for d = 1:(p-1)
    %Draw just one sample
beta = mvnrnd(mean_beta_d{d},Sigma_beta_d{d})';
indicator = binornd(1,final_w_vector_d{d});

indicatorbeta = beta .* indicator;
 %MATLAB uses the scale parameter version of gamma distribution.  To get
  %the inverse, take the reciprocal of the rate parameter,
  %because the sigma^2 ~ InvGamma and 1/(sigma^2) ~ Gamma relationship is
  %dependent on the gamma distribution having a rate parameter.

inverse_sigma_squared_d = gamrnd((capital_A + n/2), 1/s_inverseGamma_d(d));

inverse_sigma_squared_vector(d) = inverse_sigma_squared_d;

L_diag_exclude_p = sqrt(inverse_sigma_squared_d);

%Find the L using the most updated beta and inverse_sigma2
L_offdiag = -indicatorbeta.*L_diag_exclude_p;  %Don't need the pth element for the offdiagonals

L_offdiag_cell{d} = L_offdiag;

end %end of d = 1:(p-1)

%sample the pth inverse sigma_squared

last_inverse_sigma_squared = gamrnd((capital_A + n/2), 1/s_inverseGamma_p);

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

Omega_matrix = lower_triangular_matrix*lower_triangular_matrix';  %This is looking good!!

Omega_matrix_VB_samples(:,:,samps) = Omega_matrix;
end

end