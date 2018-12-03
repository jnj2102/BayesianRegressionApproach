%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=300, p=100, sparsity = twopercent

    load('CholeskyDecomp_p100_n300_twopercent_prior_FFN.mat');
    
    %Run this script to do the matrix inversion
    
  %  run('run_to_compile_and_test.m')

   c_tune_list = [.1;  1;  10];

   [num_elements, ~] = size(c_tune_list);
   
	total_time_n300_p100_twopercent = zeros(num_elements,1, reps);
    mean_Z_Bayes_est_n300_p100_twopercent = cell([num_elements,1, reps]);
    edge_matrix_n300_p100_twopercent = cell([num_elements,1, reps]);
    Sigma_Bayes_est_n300_p100_twopercent = cell([num_elements,1, reps]);
    Omega_Bayes_est_n300_p100_twopercent =  cell([num_elements,1, reps]);
        
    SP_matrix_bglasso =  zeros(num_elements,1, reps);
    SE_matrix_bglasso =  zeros(num_elements,1, reps);
    MCC_matrix_bglasso =  zeros(num_elements,1, reps);
     TP_bglasso_n300_p100_twopercent = zeros(num_elements,1, reps);
    TN_bglasso_n300_p100_twopercent = zeros(num_elements,1, reps);
    FP_bglasso_n300_p100_twopercent = zeros(num_elements,1, reps);
    FN_bglasso_n300_p100_twopercent = zeros(num_elements,1, reps);
      entropy_loss_n300_p100_twopercent = zeros(num_elements,1, reps);
    bounded_loss_n300_p100_twopercent =  zeros(num_elements,1, reps);
    Frobenius_norm_precision_n300_p100_twopercent =  zeros(num_elements,1, reps);
    Frobenius_norm_covariance_n300_p100_twopercent = zeros(num_elements,1, reps);

       
        a0 = .01;
        b0 = .01;
        capital_A = 1;
    %parfor doesn't work with the inner loop indexing.

parpool(5)


parfor (iters = 71:75)
	
   c_tune_list = [.1;  1;  10];

 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    for A_index = 1:num_elements
         rng(100 + iters,'twister'); %set the seed for each replication for reproducibility

      c_tune = c_tune_list(A_index);

      [SP_bglasso, SE_bglasso, MCC_bglasso,...
   edge_matrix_bglasso, Sigma_Bayes_est,Omega_Bayes_est,...
    total_time,mean_Z_Bayes_est, TP_bglasso, TN_bglasso,...
    FP_bglasso,FN_bglasso,Frobenius_norm_precision,entropy_loss,...
    Frobenius_norm_covariance,bounded_loss] = Bsplines_CholeskyDecomp_MCMC_horseshoe_mem(n,p,...
     x_matrix_n300_p100{iters}, omega_true, sigma_true, F_mat_cell_iters{iters},...
     g_vec_cell_iters{iters},  W_mat_cell_iters{iters}, q_vec_cell_iters{iters},...
	knot_vector_cell_iters{iters}, index_1_cell_iters{iters}, index_2_cell_iters{iters},...
    inverse_variance_prior_reduced_cell_iters{iters}, mean_prior_reduced_cell_iters{iters},...
    Z_red_cell_iters{iters}, Z_two_cell_iters{iters}, initial_value_cell_iters{iters}, capital_A, a0, b0, c_tune);
	
%save all of the data for each hyperparameter setting.
    SP_matrix_bglasso(A_index, iters) = SP_bglasso;
    SE_matrix_bglasso(A_index, iters) = SE_bglasso;
    MCC_matrix_bglasso(A_index, iters) = MCC_bglasso;
    
   
   total_time_n300_p100_twopercent(A_index, iters) = total_time; %this is the
        %total time it took to run the sampling 
 entropy_loss_n300_p100_twopercent(A_index, iters) = entropy_loss; 
   bounded_loss_n300_p100_twopercent(A_index, iters) = bounded_loss; 
   Frobenius_norm_precision_n300_p100_twopercent(A_index, iters) = Frobenius_norm_precision; 
   Frobenius_norm_covariance_n300_p100_twopercent(A_index, iters) = Frobenius_norm_covariance; 


    edge_matrix_n300_p100_twopercent{A_index, iters} = edge_matrix_bglasso;

   Omega_Bayes_est_n300_p100_twopercent{A_index, iters} = Omega_Bayes_est;
    Sigma_Bayes_est_n300_p100_twopercent{A_index, iters} = Sigma_Bayes_est;
    mean_Z_Bayes_est_n300_p100_twopercent{A_index, iters} = mean_Z_Bayes_est;
    TP_bglasso_n300_p100_twopercent(A_index, iters) = TP_bglasso;
    TN_bglasso_n300_p100_twopercent(A_index, iters) = TN_bglasso;
    FP_bglasso_n300_p100_twopercent(A_index, iters) = FP_bglasso;
    FN_bglasso_n300_p100_twopercent(A_index, iters) = FN_bglasso;

    
   end
        
end

delete(gcp('nocreate'))


save('CholeskyDecomp_p100_n300_twopercent_HS_71to75.mat', '-v7.3');


    