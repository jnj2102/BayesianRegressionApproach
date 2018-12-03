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
	
	     rng(5000,'twister'); %set the seed for each replication for reproducibility

    
[optimalJ, minK, Final_AIC,F_mat_cell,g_vec_cell,...
Fconstraint_cell,RSS_gq_cell,A_mat_cell,LSS_gq_cell,W_mat_cell,...
q_vec_cell,knot_vector_cell,index_1_cell,index_2_cell,...
inverse_variance_prior_reduced_cell,mean_prior_reduced_cell,...
Z_red_cell,Z_two_cell,initial_value] = PriorBsplines_Initialvalues_corrected(n,p, data_matrix_std, mu, tau, sigma2,c,...
N_points);


  %now run the method
  
         
        capital_A = .01;
        capital_B = .01;
		tau = 1000;
		zeta_squared = 10;
        num_VB_samples = 500; %Number of VB samples
 
 sigma_true = eye(p);
 omega_true = eye(p);
         

          
      %using identity for omega_true and sigma_true even though we don't
      %know the true omega or sigma.

[~, ~, ~, total_time, Omega_Bayes_est,~,~,...
    ~,~,~,~,~,...
    mean_Z_Bayes_est, edge_matrix_estimate] = BayesNonpar_CholeskyDecomp_VariationalBayes_plugIn_mem(n,p,...
  data_matrix_std, omega_true,sigma_true, F_mat_cell, g_vec_cell,...
    W_mat_cell, q_vec_cell,...
knot_vector_cell, index_1_cell, index_2_cell,...
inverse_variance_prior_reduced_cell,...
    mean_prior_reduced_cell,...
       Z_red_cell, Z_two_cell,...
    initial_value, capital_A,capital_B,  tau, zeta_squared,num_VB_samples);



    
  save('RealData_CholeskyDecomp_VB_Wille.mat', '-v7.3');