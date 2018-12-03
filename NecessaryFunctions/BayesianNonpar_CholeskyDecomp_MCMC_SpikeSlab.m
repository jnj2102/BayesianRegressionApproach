function [SP, SE, MCC,edge_matrix_estimate,...
   TN, TP, FN, FP, Omega_Bayes_est, total_time,mean_Z_Bayes_est, ...
   Frobenius_norm_precision,entropy_loss,...
    Frobenius_norm_covariance,bounded_loss, Sigma_Bayes_est] = BayesianNonpar_CholeskyDecomp_MCMC_SpikeSlab(n,p,...
     x_matrix, omega_true, sigma_true,F_mat_cell, g_vec_cell, W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell, inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, initial_value, capital_A, capital_B, tau, zeta_squared)

%Function to run the B-splines approach with the Cholesky decomposition for
%the Bayesian nonparanormal.

%Author: Jami Jackson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indmx = reshape(1:p^2,p,p); 
upperind = indmx(triu(indmx,1)>0);  %does not include the diagonal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Evaluate the function 

basis_full = cell([p,1]);
basis_red = cell([p,1]); 
basis_two = cell([p,1]);
  
Bayes_est_cell_initial = zeros([n,p]);
    

ind_noi_all = zeros(p-1,p);

for i = 1:p
       if i==1  
       ind_noi = (2:p)'; 
      elseif i==p
       ind_noi = (1:p-1)'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
       
end
   
%%%%%Create the basis reduced and basis two matrices to use everywhere

for d = 1:p 

%find the basis function evaluated at the gridpoints for each x

    basis_full{d} = bspline_basismatrix(4, knot_vector_cell{d}, x_matrix(:,d)); 
 
%remove the two rows of basis matrix using the two indices that match the
%thetas that we can solve for.

basis_red{d} = basis_full{d};
basis_red{d}(:, [index_1_cell{d},index_2_cell{d}]) = [];  

basis_two{d} = basis_full{d}(:, [index_1_cell{d},index_2_cell{d}]);
 
end



for d = 1:p 

%the final sample is already the dth final sample
   [Bayes_est] = estfunction(mean(initial_value{d}, 2),...  
                   basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d});
               
               Bayes_est_cell_initial(:,d)  = Bayes_est;
               
end

%initial values for the Cholesky algorithm


inverse_sigma_squared_initial_chain = ones([p,1]);

mean_initial_chain = mean(Bayes_est_cell_initial)';

Z_initial_tmp = Bayes_est_cell_initial - mean_initial_chain';

 %Tune the rho and find the initial indicators for the algorithm 
 %to sample the precision matrix using spike
%and slab prior

[final_rho_cell, indicator_vector_initial_chain] = TuneRhoandInitialW(Z_initial_tmp,zeta_squared,...
    capital_A, capital_B, tau, p,n);


%when p=n or most likely when p>n, the matrix may not be positive definite,
%so I will use the nearest positive definite matrix since these are just
%initial values anyway.    

if (p == n || p > n) == 1
    
    sigma_initial_chain = nearestSPD(cov(Bayes_est_cell_initial));
    omega_initial_chain = nearestSPD(invChol_mex(sigma_initial_chain)); %I could use the equation of an inverse of matrix 

else 
    sigma_initial_chain = cov(Bayes_est_cell_initial);
    omega_initial_chain = invChol_mex(sigma_initial_chain);
end

                                                %instead of finding the
                                                %inverse directly.
   
 
 


    burnin = 5000;
        
    
    num_samples = 15000; 

     Omega_sims_chain_tmp = zeros([p,p,num_samples]); 

    Sig_sims_chain_tmp = zeros([p,p,num_samples]);
    Bayes_est_z_iters =  zeros([n,p,num_samples]);
    

    
    thetas_samples_cell = initial_value;

    
        Omega_sparse = sparse(omega_initial_chain);

    Sig_sparse = sparse(sigma_initial_chain);
    
    mean_posterior = mean_initial_chain;
    
    Y_matrix = Bayes_est_cell_initial;
    
   indicator_vector_cell = indicator_vector_initial_chain;
    inverse_sigma_squared_vector = inverse_sigma_squared_initial_chain;
    
 
    
%Step 5 Truncated Multivariate Normal sampling using HMC
tic
for iters_num = 1:num_samples 
    rng(2000 + iters_num,'twister'); %set the seed for each replication for reproducibility

    [thetas_samples_cell,Y_matrix] = NonparanormalBsplinesPrior(F_mat_cell,thetas_samples_cell,g_vec_cell,...
    ind_noi_all,W_mat_cell, q_vec_cell, mean_posterior, Omega_sparse,...
    inverse_variance_prior_reduced_cell, mean_prior_reduced_cell, n, Z_red_cell,...
    Z_two_cell, Y_matrix,basis_red, basis_two,p);

[Z,mean_posterior] = samplemunew_Gibbs(Y_matrix, p, n, Sig_sparse);

%samplemu uses sigma initial not omega initial
Bayes_est_z_iters(:,:, iters_num) = Z;


[inverse_sigma_squared_vector, indicator_vector_cell,Omega_matrix] = CholeskyDecomp_MCMC_SpikeSlab(Z, p,n, indicator_vector_cell,...
    inverse_sigma_squared_vector, final_rho_cell, zeta_squared,...
    capital_A, capital_B);

Omega_sparse = sparse(Omega_matrix);
Omega_sparse = 1/2*(Omega_sparse + Omega_sparse'); %do this for numerical stability.
    
%take the inverse to get sigma
Sig = invChol_mex(Omega_matrix); 
Sig_sparse = sparse(Sig);
Sig_sparse = 1/2*(Sig_sparse + Sig_sparse');





Omega_sims_chain_tmp(:,:,iters_num) = Omega_sparse; 
Sig_sims_chain_tmp(:,:,iters_num) = Sig_sparse;


end %end of number of iterations   
  total_time=  toc;  %%END OF MCMC %%

%Do the burnin for the Z's.  These Z's are Z ~ N(0, Sigma). Not the edge
%matrix.
complete_Z_Bayes_est_sims = Bayes_est_z_iters(:,:,burnin+1:end);

%Do burnin for omega

complete_Omega_total_sims = Omega_sims_chain_tmp(:,:, burnin+1:end); 


complete_Sig_total_sims = Sig_sims_chain_tmp(:,:,  burnin+1:end);
   

%these are the Z ~ N(0, Sigma)

%Summaries
%Bayes estimate of Omega
Omega_Bayes_est = mean(complete_Omega_total_sims,3);

%Bayes estimate of Z transformed variables
mean_Z_Bayes_est = mean(complete_Z_Bayes_est_sims,3);


%Entropy loss
 entropy_loss = trace(Omega_Bayes_est*sigma_true) -...
     log(det(Omega_Bayes_est*sigma_true)) - p;
 
    
 
Sigma_Bayes_est = mean(complete_Sig_total_sims,3);  


 %Frobenius norm

    %Find the Frobenius Norm 
    Frobenius_norm_precision = trace((Omega_Bayes_est -...
        omega_true)'*(Omega_Bayes_est -...
        omega_true));
    
        Frobenius_norm_covariance = trace((Sigma_Bayes_est -...
        sigma_true)'*(Sigma_Bayes_est -...
        sigma_true));
    
        %Find the bounded loss
    
    bounded_loss = 1/(p^2) * sum(sum(abs(Omega_Bayes_est -...
        omega_true)));

%Omega_L1 = invChol_mex(mean(complete_Sig_total_sims,3));  %Bayes estimator using L1


%Find L1 loss

   % L1_loss = trace(Omega_L1*sigma_true) - log(det(Omega_L1*sigma_true)) - p;
    
%Now find MCC, SP, and SE


%True edge matrix for the model 
edge_matrix_true = zeros(p,p); %1 is edge is present 0 is edge is absent.

for j = 1:p
    for k = 1:p
        if omega_true(j,k) == 0
            edge_matrix_true(j,k) = 0;
        else
            edge_matrix_true(j,k) = 1;
        end
    end
end
  
%Convert the posterior sparse precision matrices to edge matrices and then
%find the edge matrix estimate using the median probability model

[~,~,num_replicates] = size(complete_Omega_total_sims);
 
complete_edge_matrix_total_sims = zeros([p,p, num_replicates]);
 
for replicates = 1:num_replicates
                  
    Omega_matrix_temp = complete_Omega_total_sims(:,:,replicates);

    edge_matrix = double(logical(Omega_matrix_temp));
 
complete_edge_matrix_total_sims(:,:,replicates) = edge_matrix;
 
    
end



%Median probability model
 edge_matrix_estimate = mean(complete_edge_matrix_total_sims,3)>.5;  %put .5 here because half the time we get 1, have the time 0

 
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication

%logical converts nonzeros to one and zeros to zero.

TP_matrix = edge_matrix_estimate(upperind) == 1  & edge_matrix_true(upperind) == 1;
TP = sum(TP_matrix(:)); %the colon sums all elements in the matrix


TN_matrix = edge_matrix_estimate(upperind) == 0 & edge_matrix_true(upperind) == 0;
TN = sum(TN_matrix(:)); 

FP_matrix = edge_matrix_estimate(upperind) == 1  & edge_matrix_true(upperind) == 0;
FP = sum(FP_matrix(:)); 

FN_matrix = edge_matrix_estimate(upperind) == 0 & edge_matrix_true(upperind) == 1;
FN = sum(FN_matrix(:)); 
%Find Specificity, Sensitivity, Matthews Correlation Coefficient 

SP = TN/(TN + FP);

SE = TP/(TP + FN);

MCC = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
  





end

