%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=300, p=100, sparsity = circle

    load('CholeskyDecomp_p100_n300_circle_prior_FFN.mat');
    
    %Run this script to do the matrix inversion
    
   % run('run_to_compile_and_test.m')

  num_elements = 1;
   
	total_time_n300_p100_circle = zeros(num_elements,1, reps);
    mean_Z_Bayes_est_n300_p100_circle = cell([num_elements,1, reps]);
    edge_matrix_n300_p100_circle = cell([num_elements,1, reps]);
    Sigma_Bayes_est_n300_p100_circle = cell([num_elements,1, reps]);
    Omega_Bayes_est_n300_p100_circle =  cell([num_elements,1, reps]);
        
    SP_matrix =  zeros(num_elements,1, reps);
    SE_matrix =  zeros(num_elements,1, reps);
    MCC_matrix =  zeros(num_elements,1, reps);
    entropy_loss_n300_p100_circle = zeros(num_elements,1, reps);
    bounded_loss_n300_p100_circle =  zeros(num_elements,1, reps);
    Frobenius_norm_precision_n300_p100_circle =  zeros(num_elements,1, reps);
    Frobenius_norm_covariance_n300_p100_circle = zeros(num_elements,1, reps);

     TP_n300_p100_circle = zeros(num_elements,1, reps);
    TN_n300_p100_circle = zeros(num_elements,1, reps);
    FP_n300_p100_circle = zeros(num_elements,1, reps);
    FN_n300_p100_circle = zeros(num_elements,1, reps);
       
	   capital_A = .01;
        capital_B = .01;
		tau = 1000;
		zeta_squared = 10;
		
%parpool(5)

for iters = 3:3
	
    rng(100 + iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
[SP, SE, MCC,edge_matrix_estimate,...
   TN, TP, FN, FP, Omega_Bayes_est, total_time,mean_Z_Bayes_est, ...
   Frobenius_norm_precision,entropy_loss,...
    Frobenius_norm_covariance,bounded_loss, Sigma_Bayes_est] = BayesianNonpar_CholeskyDecomp_MCMC_SpikeSlab_mem(n,p,...
     x_matrix_n300_p100{iters}, omega_true, sigma_true, F_mat_cell_iters{iters}, g_vec_cell_iters{iters},...
     W_mat_cell_iters{iters}, q_vec_cell_iters{iters},...
knot_vector_cell_iters{iters}, index_1_cell_iters{iters}, index_2_cell_iters{iters}, inverse_variance_prior_reduced_cell_iters{iters},...
 mean_prior_reduced_cell_iters{iters},...
Z_red_cell_iters{iters}, Z_two_cell_iters{iters},...
initial_value_cell_iters{iters}, capital_A, capital_B, tau, zeta_squared);
          

	
%save all of the data for each hyperparameter setting.
    SP_matrix( iters) = SP;
    SE_matrix( iters) = SE;
    MCC_matrix( iters) = MCC;
    
   
   total_time_n300_p100_circle( iters) = total_time; %this is the
        %total time it took to run the sampling 
entropy_loss_n300_p100_circle(iters) = entropy_loss; 
   bounded_loss_n300_p100_circle(iters) = bounded_loss; 
   Frobenius_norm_precision_n300_p100_circle(iters) = Frobenius_norm_precision; 
   Frobenius_norm_covariance_n300_p100_circle(iters) = Frobenius_norm_covariance; 


    edge_matrix_n300_p100_circle{ iters} = edge_matrix_estimate;

   Omega_Bayes_est_n300_p100_circle{ iters} = Omega_Bayes_est;
    Sigma_Bayes_est_n300_p100_circle{ iters} = Sigma_Bayes_est;
    mean_Z_Bayes_est_n300_p100_circle{ iters} = mean_Z_Bayes_est;
    TP_n300_p100_circle( iters) = TP;
    TN_n300_p100_circle( iters) = TN;
    FP_n300_p100_circle( iters) = FP;
    FN_n300_p100_circle( iters) = FN;

    

        
end

%delete(gcp('nocreate'))

save('CholeskyDecomp_p100_n300_circle_SS_3to3.mat', '-v7.3');


    