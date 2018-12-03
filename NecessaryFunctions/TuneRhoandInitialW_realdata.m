function [final_rho_cell, final_w_vector_cell] = TuneRhoandInitialW_realdata(Z,zeta_squared,...
    capital_A, capital_B, tau, p,n)

%Tune the rho and find the initial values for w for the VB algorithm for
%the Cholesky decomposition
%Author: Jami Jackson Mulgrave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_value = 100;
P_value = 50;
exp_value = exp(-0.5*sqrt(n));
lambda_j = linspace(-15,5,50);
c_j = 1/p*linspace(0.1,10,50);
exp_lambda_j = exp(lambda_j)./(1+exp(lambda_j));

w_vector_cell = cell([p-1,1]);
rho_cell = cell([p-1,1]);

%Initialize w vector with zeros
for m = 1:(p-1) 
    
    w_temp = zeros([p-m,1]);
    w_vector_cell{m} = w_temp;
    row_offdiag = (m+1:p)';
    rho_j_temp = (exp_lambda_j./sqrt(row_offdiag)) .* c_j; %the rho gets updated
    rho_cell{m} = rho_j_temp;

end


for d = 1:(p-1)
    L_value = -inf; %has to restart for each regression

row_offdiag = (d+1:p)';  %Need the loop indexed
[dim_w_vector, ~] = size(row_offdiag);

        %rho/sqrt(row_index)
        rho = (exp_value/(1+exp_value))./(sqrt(row_offdiag)) .* 1/p;

        Z_designMat = Z(:, (d+1):end);
        Z_outcome = Z(:,d);

        lower_bound_temp = zeros([dim_w_vector, 1]);

  for i = 1:max(p, P_value)
      

        %I need to call this w_vector_cell back for each i because I'm
        %adding the j each time.
        w_vector = w_vector_cell{d}; %this updates after each L update

        %the w_vector is only relevant for the p-1 regressors
        for j = 1:dim_w_vector
        w_vector_j = w_vector;
        w_vector_j(j) = 1;
        

   lower_bound_temp(j) = LowerboundMarginalLikelihoodCholeskyDecomp(Z_designMat,...
    Z_outcome, zeta_squared, capital_A, capital_B, tau, w_vector_j, n, rho);
    
        end
        
      
        %find the k which is the maximum for the dth regression (for
        %p-1). We don't need the pth regression because it isn't affected
        %by the w vector - the w vector is constant with respect to the pth
        %regression

[max_lower_bound_k, k_index] = max(lower_bound_temp);


if max_lower_bound_k > L_value
    
    L_value = max_lower_bound_k;
    w_vector(k_index) = 1;
    
    w_vector_cell{d} = w_vector;
    continue;
else 
    break;

end
    
  end
  
  
%Now tune the rho
lower_bound_temp = zeros([P_value, 1]);
%L_value = -inf*ones([p-1,1]); 


%L_value = -inf; %has to restart for each regression

for m = 1:M_value
    
   rho_j_temp = rho_cell{d};
    
    for j = 1:P_value
        
        rho_j = rho_j_temp(:,j);
        
        lower_bound_temp(j) = LowerboundMarginalLikelihoodCholeskyDecomp(Z_designMat,...
    Z_outcome, zeta_squared, capital_A, capital_B, tau, w_vector,...
    n, rho_j);

        
    end


        %find the k which is the maximum for the dth regression (for
        %p-1). We don't need the pth regression because it isn't affected
        %by the rho vector - the rho vector is constant with respect to the pth
        %regression

[max_lower_bound_k, k_index] = max(lower_bound_temp);


if max_lower_bound_k >= L_value
    
    L_value = max_lower_bound_k;
    final_rho = rho_j_temp(:,k_index); %select the final rho
    
end

  %  L_value = -inf; %has to restart for each regression

for j = 1:dim_w_vector
     w_vector_j_zero = w_vector;
      w_vector_j_zero(j) = 0;
      
    lower_bound_zero = LowerboundMarginalLikelihoodCholeskyDecomp(Z_designMat,...
    Z_outcome, zeta_squared, capital_A, capital_B, tau, w_vector_j_zero,...
    n, final_rho );

     w_vector_j_one = w_vector;
      w_vector_j_one(j) = 1;
      
    lower_bound_one = LowerboundMarginalLikelihoodCholeskyDecomp(Z_designMat,...
    Z_outcome, zeta_squared, capital_A, capital_B, tau, w_vector_j_one,...
    n, final_rho);

    w_vector_j_both = [w_vector_j_zero(j),w_vector_j_one(j)];

    [max_bound, k_index_max] = max([lower_bound_zero, lower_bound_one]);
    
    if max_bound > L_value
        L_value = max_bound;
       w_vector(j) = w_vector_j_both(k_index_max);
    
    end

end

if max_bound <= L_value
    final_rho_cell{d} = final_rho;
    final_w_vector_cell{d} = w_vector;
    break;
    
end

end
    

end


    





end








