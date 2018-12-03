%%Real data application for Bayesian Regression paper
% Estimating transformation functions.
%Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear;

load('Bsplines_paper_realdata_initialdata.mat');

%log transform the data and standardize the each gene

data_matrix_log = log(data_matrix)';

data_matrix_std = (data_matrix_log - min(data_matrix_log))./(max(data_matrix_log) - min(data_matrix_log));


%Now try my method out - data_matrix_std is my "x_matrix"

[n,p] = size(data_matrix_std);

%these (mu, tau, sigma2) are for my prior for the transformation functions
mu = 1; %I'm just giving it a mean of 1 for now but it can be any constant
tau = 1; %I'm giving it a sd of 1 for now but it can be any constant

sigma2 = 1; %I'm making the variance 1 for now but it can be any constant

N_points = 20;

c = [0;1]; %vector of linear constraints
 %First find the optimal J for the replication, then the prior mean and
    %variance and finally the initial values
    
	     rng(1000,'twister'); %set the seed for each replication for reproducibility

[optimalJ, minK, Final_AIC,F_mat_cell,g_vec_cell,...
Fconstraint_cell,RHS_gq_cell,A_mat_cell,LHS_gq_cell,W_mat_cell,...
q_vec_cell,knot_vector_cell,index_1_cell,index_2_cell,...
inverse_variance_prior_reduced_cell,mean_prior_reduced_cell,...
Z_red_cell,Z_two_cell,initial_value] = PriorBsplines_Initialvalues_corrected(n,p, data_matrix_std, mu, tau, sigma2,c,...
N_points);


  %now run the method
  
   c_tune_list = [.1;  1;  10];

   [num_elements, ~] = size(c_tune_list);
   
	total_time_realdata = zeros(num_elements,1, 1);
    mean_Z_Bayes_est_realdata = cell([num_elements,1, 1]);
    edge_matrix_realdata = cell([num_elements,1, 1]);
    Sigma_Bayes_est_realdata = cell([num_elements,1, 1]);
    Omega_Bayes_est_realdata =  cell([num_elements,1, 1]);
        
       
        a0 = .01;
        b0 = .01;
        capital_A = 1;

		 
 sigma_true = eye(p);
 omega_true = eye(p);
         
    for A_index = 1:num_elements
     rng(5000,'twister'); %set the seed for each replication for reproducibility

      c_tune = c_tune_list(A_index);
          
      %using identity for omega_true and sigma_true even though we don't
      %know the true omega or sigma.

    [~, ~, ~,...
   edge_matrix_bglasso, Sigma_Bayes_est,Omega_Bayes_est,...
    total_time,mean_Z_Bayes_est, ~, ~,...
    ~,~,~,~,...
    ~,~] = Bsplines_CholeskyDecomp_MCMC_horseshoe_mem(n,p,...
     data_matrix_std, omega_true, sigma_true, F_mat_cell,...
     g_vec_cell,  W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell,...
    inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, initial_value, capital_A, a0, b0, c_tune);


%save all of the data for each hyperparameter setting.

   
   total_time_realdata(A_index) = total_time; %this is the

        %total time it took to run the sampling 
    edge_matrix_realdata{A_index} = edge_matrix_bglasso;

   Omega_Bayes_est_realdata{A_index} = Omega_Bayes_est;
    Sigma_Bayes_est_realdata{A_index} = Sigma_Bayes_est;
    mean_Z_Bayes_est_realdata{A_index} = mean_Z_Bayes_est;

    
    end
        
   
     
     for A_index = 1:num_elements
               
        Omega_bayes_est_temp =  Omega_Bayes_est_realdata{A_index};
         edge_matrix_temp = edge_matrix_realdata{A_index};
        
       mean_Z_Bayes_est_temp = mean_Z_Bayes_est_realdata{A_index};
       
       %%%Find the MLE using a convex optimization problem for eBIC

     %%% MLE with Bglasso
     
[BIC]= FindBIC(Omega_bayes_est_temp,edge_matrix_temp,mean_Z_Bayes_est_temp,n,p);

     BIC_matrix(A_index) = BIC;
   
     end
    
     [minBIC, ~ ] = min(min(BIC_matrix(:)));

[rowMinBIC, ~] = find(BIC_matrix(:) == minBIC);

%save the data with the minimum BIC for further analysis
 
BIC_matrix_finalanalysis =   BIC_matrix(rowMinBIC);

 total_time_finalanalysis =   total_time_realdata(rowMinBIC);
 Omega_Bayes_est_finalanalysis =   Omega_Bayes_est_realdata{rowMinBIC};
Sigma_Bayes_est_finalanalysis = Sigma_Bayes_est_realdata{rowMinBIC};
edge_matrix_finalanalysis =  edge_matrix_realdata{rowMinBIC};

   
  mean_Z_Bayes_est_finalanalysis = mean_Z_Bayes_est_realdata{rowMinBIC};
  
save('RealData_CholeskyDecomp_HS_Wille.mat', '-v7.3');