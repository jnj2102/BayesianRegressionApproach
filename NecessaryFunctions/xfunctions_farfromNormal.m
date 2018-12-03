function [x_matrix, MLE_ald,MLE_ev, MLE_stable] = xfunctions_farfromNormal(p, y_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate the x matrix  using MLE
%Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the maximum likelihood estimate for each column of Y
%(corresponding to the predictor) using the distribution that corresponds
%to the cdf

MLE_ald =  cell([p,1]);

MLE_ev = cell([p,1]);
MLE_stable = cell([p,1]);

x_matrix = [];

d = 0;

while (1)    
 d = d+1;
    
        if d > p
            break
        else

            %Use the ald package in R
            
csvwrite('y_true_d.csv',y_true(:, d));
system('source /usr/local/apps/R/R-312.csh; R CMD BATCH convert_to_ald.R') %R must be in the path for this to run
ald_cdf_xmat = csvread('ald_CDF_for_xmatrix.csv',1,1);
MLE_ald{d} = csvread('ald_MLE.csv',1,1); %mu sigma and skewness
            
x_matrix = horzcat(x_matrix, ald_cdf_xmat);

        end
        
  d = d+1;
  
  if d > p
       break
  else
         
MLE_ev{d} = mle(y_true(:, d),'distribution','Extreme Value');
x_matrix = horzcat(x_matrix, cdf('Extreme Value', y_true(:, d), MLE_ev{d}(1), MLE_ev{d}(2)));
  end
    
  d = d+1;
  
  if d > p
      break
  else
%the beta is on the boundary - this may or may not be ok for stable.
MLE_stable{d} = mle(y_true(:, d),'distribution', 'stable');
x_matrix = horzcat(x_matrix, cdf('stable', y_true(:, d), MLE_stable{d}(1),...
    MLE_stable{d}(2), MLE_stable{d}(3), MLE_stable{d}(4)));
  end

   
end








end