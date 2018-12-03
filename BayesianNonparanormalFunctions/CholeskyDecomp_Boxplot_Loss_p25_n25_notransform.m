%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the CholeskyDecomp paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AR2 p25 n25 SS 

load('CholeskyDecomp_p25_n25_AR2_SS_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_AR2, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_AR2, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_AR2, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_AR2_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [table_p25_n25_AR2_SS];


clearvars -except combine_tables

%AR2 p25 n25 VB 

load('CholeskyDecomp_p25_n25_AR2_VB_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_AR2, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_AR2, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_AR2, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_AR2_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_AR2_VB];


clearvars -except combine_tables

%AR2 p25 n25 HS 

load('CholeskyDecomp_p25_n25_AR2_HS_notransform_final.mat');

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

table_p25_n25_AR2_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_AR2_HS];


clearvars -except combine_tables

%tenpercent p25 n25 SS 

load('CholeskyDecomp_p25_n25_tenpercent_SS_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_tenpercent, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_tenpercent, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_tenpercent, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_tenpercent_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_tenpercent_SS];


clearvars -except combine_tables

%tenpercent p25 n25 VB 

load('CholeskyDecomp_p25_n25_tenpercent_VB_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_tenpercent, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_tenpercent, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_tenpercent, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_tenpercent_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_tenpercent_VB];


clearvars -except combine_tables

%tenpercent p25 n25 HS 

load('CholeskyDecomp_p25_n25_tenpercent_HS_notransform_final.mat');

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

table_p25_n25_tenpercent_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_tenpercent_HS];


clearvars -except combine_tables

%circle p25 n25 SS 

load('CholeskyDecomp_p25_n25_circle_SS_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_circle, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_circle, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_circle, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_circle_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p25_n25_circle_SS];


clearvars -except combine_tables

%circle p25 n25 VB 

load('CholeskyDecomp_p25_n25_circle_VB_notransform.mat');

ELoss_value = reshape(entropy_loss_n25_p25_circle, [reps,1]);
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = reshape(bounded_loss_n25_p25_circle, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = reshape(Frobenius_norm_precision_n25_p25_circle, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n25_circle_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_circle_VB];


clearvars -except combine_tables

%circle p25 n25 HS 

load('CholeskyDecomp_p25_n25_circle_HS_notransform_final.mat');

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

table_p25_n25_circle_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n25_circle_HS];


clearvars -except combine_tables

%Try writing table as a csv file to read into R

writetable(combine_tables,'CholeskyDecomp_Boxplot_Loss_p25_n25_notransform.csv') 

