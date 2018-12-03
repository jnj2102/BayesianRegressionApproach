function [SP, SE, MCC,edge_matrix_estimate,...
   TN, TP, FN, FP, Omega_Bayes_est, total_time,mean_Z_Bayes_est, ...
   Frobenius_norm_precision,entropy_loss,...
    Frobenius_norm_covariance,bounded_loss, Sigma_Bayes_est] = NoTransform_CholeskyDecomp_MCMC_SpikeSlab_mem(n,p,...
     x_matrix, omega_true, sigma_true, capital_A, capital_B, tau, zeta_squared)

%Function to run the B-splines approach with the Cholesky decomposition for
%the Bayesian nonparanormal.

%Author: Jami Jackson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indmx = reshape(1:p^2,p,p); 
upperind = indmx(triu(indmx,1)>0);  %does not include the diagonal


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
   
%initial values for the Cholesky algorithm


inverse_sigma_squared_initial_chain = ones([p,1]);

mean_initial_chain = mean(x_matrix)';

Z_initial_tmp = x_matrix - mean_initial_chain';

 %Tune the rho and find the initial indicators for the algorithm 
 %to sample the precision matrix using spike
%and slab prior

[final_rho_cell, indicator_vector_initial_chain] = TuneRhoandInitialW(Z_initial_tmp,zeta_squared,...
    capital_A, capital_B, tau, p,n);


%when p=n or most likely when p>n, the matrix may not be positive definite,
%so I will use the nearest positive definite matrix since these are just
%initial values anyway.    

if (p == n || p > n) == 1
    
    sigma_initial_chain = nearestSPD(cov(x_matrix));

else 
    sigma_initial_chain = cov(x_matrix);
end

                                                %instead of finding the
                                                %inverse directly.
   
 
 


    burnin = 5000;
        
    
    num_samples = 15000; 



    complete_Z_Bayes_est_sims = zeros([n,p,num_samples - burnin]);
    complete_Omega_total_sims = zeros([p,p,num_samples - burnin]);
     complete_Sig_total_sims = zeros([p,p,num_samples - burnin]);
    
    
    Sig_sparse = sparse(sigma_initial_chain);
            
   indicator_vector_cell = indicator_vector_initial_chain;
    inverse_sigma_squared_vector = inverse_sigma_squared_initial_chain;
    
 
    
%Step 5 Truncated Multivariate Normal sampling using HMC
tic
for iters_num = 1:num_samples 
    rng(2000 + iters_num,'twister'); %set the seed for each replication for reproducibility

% %update my thetas with all of the last drawn thetas.

[Z,~] = samplemunew_Gibbs(x_matrix, p, n, Sig_sparse);

%samplemu uses sigma initial not omega initial


[inverse_sigma_squared_vector, indicator_vector_cell,Omega_matrix] = CholeskyDecomp_MCMC_SpikeSlab(Z, p,n, indicator_vector_cell,...
    inverse_sigma_squared_vector, final_rho_cell, zeta_squared,...
    capital_A, capital_B);

Omega_sparse = sparse(Omega_matrix);
Omega_sparse = 1/2*(Omega_sparse + Omega_sparse'); %do this for numerical stability.
    
%take the inverse to get sigma
Sig = invChol_mex(Omega_matrix); 
Sig_sparse = sparse(Sig);
Sig_sparse = 1/2*(Sig_sparse + Sig_sparse');



%Only save the samples after burnin
if iters_num > burnin
    %create the corresponding index
    
   save_index = iters_num - burnin;

complete_Sig_total_sims(:,:,save_index) = Sig_sparse;


complete_Omega_total_sims(:,:,save_index) = Omega_sparse;

complete_Z_Bayes_est_sims(:,:,save_index) = Z; %these cost some time

end



end %end of number of iterations   
  total_time=  toc;  %%END OF MCMC %%

%Do the burnin for the Z's.  These Z's are Z ~ N(0, Sigma). Not the edge
%matrix.
% complete_Z_Bayes_est_sims = Bayes_est_z_iters(:,:,burnin+1:end);
% 
% %Do burnin for omega
% 
% complete_Omega_total_sims = Omega_sims_chain_tmp(:,:, burnin+1:end); 
% 
% 
% complete_Sig_total_sims = Sig_sims_chain_tmp(:,:,  burnin+1:end);
%    
% 
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

