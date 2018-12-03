%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the CholeskyDecomp paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AR2 p50 n100 SS 

 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p50_n100_AR2_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n100_p50_AR2_tmp = entropy_loss_n100_p50_AR2(:);
bounded_loss_n100_p50_AR2_tmp = bounded_loss_n100_p50_AR2(:);
Frobenius_norm_precision_n100_p50_AR2_tmp = Frobenius_norm_precision_n100_p50_AR2(:);

entropy_loss_n100_p50_AR2 = entropy_loss_n100_p50_AR2_tmp(ssIters:ssIters2);
bounded_loss_n100_p50_AR2 = bounded_loss_n100_p50_AR2_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n100_p50_AR2 = Frobenius_norm_precision_n100_p50_AR2_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n100_p50_AR2];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n100_p50_AR2];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n100_p50_AR2];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 


ELoss_value = reshape(ELoss_value_temp, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(BLoss_value_temp, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(FrobLoss_value_temp, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);


Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR2_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [table_p50_n100_AR2_SS];


clearvars -except combine_tables

%AR2 p50 n100 VB 

load('CholeskyDecomp_p50_n100_AR2_VB.mat');

ELoss_value = reshape(entropy_loss_n100_p50_AR2, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n100_p50_AR2, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n100_p50_AR2, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR2_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR2_VB];


clearvars -except combine_tables

%AR2 p50 n100 HS 

load('CholeskyDecomp_p50_n100_AR2_HS_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR2_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR2_HS];


clearvars -except combine_tables

%fivepercent p50 n100 SS 

 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p50_n100_fivepercent_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n100_p50_fivepercent_tmp = entropy_loss_n100_p50_fivepercent(:);
bounded_loss_n100_p50_fivepercent_tmp = bounded_loss_n100_p50_fivepercent(:);
Frobenius_norm_precision_n100_p50_fivepercent_tmp = Frobenius_norm_precision_n100_p50_fivepercent(:);

entropy_loss_n100_p50_fivepercent = entropy_loss_n100_p50_fivepercent_tmp(ssIters:ssIters2);
bounded_loss_n100_p50_fivepercent = bounded_loss_n100_p50_fivepercent_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n100_p50_fivepercent = Frobenius_norm_precision_n100_p50_fivepercent_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n100_p50_fivepercent];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n100_p50_fivepercent];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n100_p50_fivepercent];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 


ELoss_value = reshape(ELoss_value_temp, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(BLoss_value_temp, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(FrobLoss_value_temp, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_fivepercent_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p50_n100_fivepercent_SS];


clearvars -except combine_tables

%fivepercent p50 n100 VB 

load('CholeskyDecomp_p50_n100_fivepercent_VB.mat');

ELoss_value = reshape(entropy_loss_n100_p50_fivepercent, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n100_p50_fivepercent, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n100_p50_fivepercent, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_fivepercent_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_fivepercent_VB];


clearvars -except combine_tables

%fivepercent p50 n100 HS 

load('CholeskyDecomp_p50_n100_fivepercent_HS_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_fivepercent_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_fivepercent_HS];


clearvars -except combine_tables

%circle p50 n100 SS 

 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p50_n100_circle_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n100_p50_circle_tmp = entropy_loss_n100_p50_circle(:);
bounded_loss_n100_p50_circle_tmp = bounded_loss_n100_p50_circle(:);
Frobenius_norm_precision_n100_p50_circle_tmp = Frobenius_norm_precision_n100_p50_circle(:);

entropy_loss_n100_p50_circle = entropy_loss_n100_p50_circle_tmp(ssIters:ssIters2);
bounded_loss_n100_p50_circle = bounded_loss_n100_p50_circle_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n100_p50_circle = Frobenius_norm_precision_n100_p50_circle_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n100_p50_circle];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n100_p50_circle];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n100_p50_circle];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 


ELoss_value = reshape(ELoss_value_temp, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(BLoss_value_temp, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(FrobLoss_value_temp, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);


Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_circle_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p50_n100_circle_SS];


clearvars -except combine_tables

%circle p50 n100 VB 

load('CholeskyDecomp_p50_n100_circle_VB.mat');

ELoss_value = reshape(entropy_loss_n100_p50_circle, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n100_p50_circle, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n100_p50_circle, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_circle_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_circle_VB];


clearvars -except combine_tables

%circle p50 n100 HS 

load('CholeskyDecomp_p50_n100_circle_HS_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_circle_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_circle_HS];


clearvars -except combine_tables

%Try writing table as a csv file to read into R

writetable(combine_tables,'CholeskyDecomp_Boxplot_Loss_p50_n100.csv') 

