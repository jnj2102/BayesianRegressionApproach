%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=300, p=100, sparsity = AR2

    load('CholeskyDecomp_p100_n300_AR2_prior_FFN.mat');
    
    %Run this script to do the matrix inversion
    
   % run('run_to_compile_and_test.m')

  num_elements = 1;
   
	total_time_n300_p100_AR2 = zeros(num_elements,1, reps);
    mean_Z_Bayes_est_n300_p100_AR2 = cell([num_elements,1, reps]);
    edge_matrix_n300_p100_AR2 = cell([num_elements,1, reps]);
    Sigma_Bayes_est_n300_p100_AR2 = cell([num_elements,1, reps]);
    Omega_Bayes_est_n300_p100_AR2 =  cell([num_elements,1, reps]);
        
    SP_matrix =  zeros(num_elements,1, reps);
    SE_matrix =  zeros(num_elements,1, reps);
    MCC_matrix =  zeros(num_elements,1, reps);
    entropy_loss_n300_p100_AR2 = zeros(num_elements,1, reps);
    bounded_loss_n300_p100_AR2 =  zeros(num_elements,1, reps);
    Frobenius_norm_precision_n300_p100_AR2 =  zeros(num_elements,1, reps);
    Frobenius_norm_covariance_n300_p100_AR2 = zeros(num_elements,1, reps);

     TP_n300_p100_AR2 = zeros(num_elements,1, reps);
    TN_n300_p100_AR2 = zeros(num_elements,1, reps);
    FP_n300_p100_AR2 = zeros(num_elements,1, reps);
    FN_n300_p100_AR2 = zeros(num_elements,1, reps);
       
	   capital_A = .01;
        capital_B = .01;
		tau = 1000;
		zeta_squared = 10;
		
parpool(8)

parfor (iters = 26:50) 
	
    rng(100 + iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    [SP, SE, MCC,edge_matrix_estimate,...
   TN, TP, FN, FP, Omega_Bayes_est, total_time,mean_Z_Bayes_est, ...
   Frobenius_norm_precision,entropy_loss,...
    Frobenius_norm_covariance,bounded_loss, Sigma_Bayes_est] = NoTransform_CholeskyDecomp_MCMC_SpikeSlab_mem(n,p,...
      x_matrix_n300_p100{iters}, omega_true, sigma_true, capital_A, capital_B, tau, zeta_squared);      

	
%save all of the data for each hyperparameter setting.
    SP_matrix( iters) = SP;
    SE_matrix( iters) = SE;
    MCC_matrix( iters) = MCC;
    
   
   total_time_n300_p100_AR2( iters) = total_time; %this is the
        %total time it took to run the sampling 
entropy_loss_n300_p100_AR2(iters) = entropy_loss; 
   bounded_loss_n300_p100_AR2(iters) = bounded_loss; 
   Frobenius_norm_precision_n300_p100_AR2(iters) = Frobenius_norm_precision; 
   Frobenius_norm_covariance_n300_p100_AR2(iters) = Frobenius_norm_covariance; 


    edge_matrix_n300_p100_AR2{ iters} = edge_matrix_estimate;

   Omega_Bayes_est_n300_p100_AR2{ iters} = Omega_Bayes_est;
    Sigma_Bayes_est_n300_p100_AR2{ iters} = Sigma_Bayes_est;
    mean_Z_Bayes_est_n300_p100_AR2{ iters} = mean_Z_Bayes_est;
    TP_n300_p100_AR2( iters) = TP;
    TN_n300_p100_AR2( iters) = TN;
    FP_n300_p100_AR2( iters) = FP;
    FN_n300_p100_AR2( iters) = FN;

    

        
end

delete(gcp('nocreate'))

save('CholeskyDecomp_p100_n300_AR2_SS_notransform_26to50.mat', '-v7.3');


    