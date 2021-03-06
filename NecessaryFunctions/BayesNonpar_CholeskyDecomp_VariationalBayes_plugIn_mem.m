function[SP, SE, MCC, total_time, Omega_Bayes_est,Frobenius_norm_precision,entropy_loss,...
    bounded_loss,TP,TN,FP,FN,...
    mean_Z_Bayes_est, edge_matrix_estimate] = BayesNonpar_CholeskyDecomp_VariationalBayes_plugIn_mem(n,p,...
     x_matrix, omega_true, sigma_true,F_mat_cell, g_vec_cell, W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell, inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, initial_value, capital_A,capital_B,  tau, zeta_squared,num_VB_samples)
%Function to do the transformation via B-spline prior and the variational
%Bayes algorithm to estimate the precision matrix.
%
%Author: Jami Jackson Mulgrave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

   [Bayes_est] = estfunction(mean(initial_value{d}, 2),...  
                   basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d});
               
               Bayes_est_cell_initial(:,d)  = Bayes_est;
               
end


%initial values for the Cholesky algorithm

mean_initial_chain = mean(Bayes_est_cell_initial)';

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
                                      
 
 

burnin = 5000;
        

num_samples = 15000; 

  complete_Z_Bayes_est_sims = zeros([n,p,num_samples - burnin]);
 
    thetas_samples_cell = initial_value;
    
    mean_posterior = mean_initial_chain;
    
    Y_matrix = Bayes_est_cell_initial;
  
  
    
        Omega_sparse = sparse(omega_initial_chain);

    Sig_sparse = sparse(sigma_initial_chain);

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

     
%Only save the samples after burnin
if iters_num > burnin
    %create the corresponding index
    
   save_index = iters_num - burnin;


complete_Z_Bayes_est_sims(:,:,save_index) = Z; %these cost some time

end



end %end of number of iterations
%End of MCMC Algorithm.

%Do the burnin for the Z's.  These Z's are Z ~ N(0, Sigma). Not the edge
%matrix.
% complete_Z_Bayes_est_sims = Bayes_est_z_iters(:,:,burnin+1:end);

%Bayes estimate of Z transformed variables
mean_Z_Bayes_est = mean(complete_Z_Bayes_est_sims,3);

 %Finding the average of the Z's is the same as finding 
 %Z = \theta_bar*B(X)- mu_bar
 
%Initialize the VB algorithm with the tuned rho and initial w for all iterations

[final_rho_cell, final_w_vector_cell] = TuneRhoandInitialW(mean_Z_Bayes_est,...
    zeta_squared,capital_A, capital_B, tau, p,n);

%Do the Variational Bayes algorithm


[Omega_matrix_VB_samples] = VBandOmega(mean_Z_Bayes_est, zeta_squared,...
    capital_A, capital_B, tau, final_w_vector_cell, n, p,final_rho_cell, num_VB_samples);
total_time = toc;
%Find the total time of MCMC + VB


%Summaries
%Bayes estimate of Omega
Omega_Bayes_est = mean(Omega_matrix_VB_samples,3);

%Entropy loss
 entropy_loss = trace(Omega_Bayes_est*sigma_true) -...
     log(det(Omega_Bayes_est*sigma_true)) - p;
 
 %Find the Frobenius Norm 
    Frobenius_norm_precision = trace((Omega_Bayes_est -...
        omega_true)'*(Omega_Bayes_est -...
        omega_true));
    
         %Find the bounded loss
    
    bounded_loss = 1/(p^2) * sum(sum(abs(Omega_Bayes_est -...
        omega_true)));   

%Find L1 loss

  %  L1_loss = trace(Omega_matrix_estimate*sigma_true) - log(det(Omega_matrix_estimate*sigma_true)) - p;

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


[~,~,num_replicates] = size(Omega_matrix_VB_samples);
 
complete_edge_matrix_total_sims = zeros([p,p, num_replicates]);
 
for replicates = 1:num_replicates
                  
    Omega_matrix_temp = Omega_matrix_VB_samples(:,:,replicates);

    edge_matrix = double(logical(Omega_matrix_temp));
 
complete_edge_matrix_total_sims(:,:,replicates) = edge_matrix;
 
    
end


%Median probability model
 edge_matrix_estimate = mean(complete_edge_matrix_total_sims,3)>.5;  %put .5 here because half the time we get 1, have the time 0

  
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication

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