%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the CholeskyDecomp paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AR2 p100 n300 SS 

 ssIters = 1;
 ssIters2 = ssIters ;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 25

load(sprintf('CholeskyDecomp_p100_n300_AR2_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n300_p100_AR2_tmp = entropy_loss_n300_p100_AR2(:);
bounded_loss_n300_p100_AR2_tmp = bounded_loss_n300_p100_AR2(:);
Frobenius_norm_precision_n300_p100_AR2_tmp = Frobenius_norm_precision_n300_p100_AR2(:);

entropy_loss_n300_p100_AR2 = entropy_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_AR2 = bounded_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_AR2 = Frobenius_norm_precision_n300_p100_AR2_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_AR2];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_AR2];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_AR2];

ssIters = ssIters2 +1;
ssIters2 = ssIters ;

 end
 

 clear ssIters ssIters2 

 ssIters = 26;
 ssIters2 = ssIters + 4;

 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_AR2_SS_%dto%d.mat', ssIters, ssIters2));

 
entropy_loss_n300_p100_AR2_tmp = entropy_loss_n300_p100_AR2(:);
bounded_loss_n300_p100_AR2_tmp = bounded_loss_n300_p100_AR2(:);
Frobenius_norm_precision_n300_p100_AR2_tmp = Frobenius_norm_precision_n300_p100_AR2(:);

entropy_loss_n300_p100_AR2 = entropy_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_AR2 = bounded_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_AR2 = Frobenius_norm_precision_n300_p100_AR2_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_AR2];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_AR2];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_AR2];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 4;

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

table_p100_n300_AR2_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [table_p100_n300_AR2_SS];


clearvars -except combine_tables

%AR2 p100 n300 VB 


 ssIters = 1;
 ssIters2 = ssIters + 24;

  ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_AR2_VB_%dto%d.mat', ssIters, ssIters2));
 
entropy_loss_n300_p100_AR2_tmp = entropy_loss_n300_p100_AR2(:);
bounded_loss_n300_p100_AR2_tmp = bounded_loss_n300_p100_AR2(:);
Frobenius_norm_precision_n300_p100_AR2_tmp = Frobenius_norm_precision_n300_p100_AR2(:);

entropy_loss_n300_p100_AR2 = entropy_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_AR2 = bounded_loss_n300_p100_AR2_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_AR2 = Frobenius_norm_precision_n300_p100_AR2_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_AR2];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_AR2];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_AR2];

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
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_AR2_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_AR2_VB];


clearvars -except combine_tables

%AR2 p100 n300 HS 

load('CholeskyDecomp_p100_n300_AR2_HS_final.mat');

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

table_p100_n300_AR2_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_AR2_HS];


clearvars -except combine_tables
%twopercent p100 n300 SS 

 ssIters = 1;
 ssIters2 = ssIters ;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 25

load(sprintf('CholeskyDecomp_p100_n300_twopercent_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n300_p100_twopercent_tmp = entropy_loss_n300_p100_twopercent(:);
bounded_loss_n300_p100_twopercent_tmp = bounded_loss_n300_p100_twopercent(:);
Frobenius_norm_precision_n300_p100_twopercent_tmp = Frobenius_norm_precision_n300_p100_twopercent(:);

entropy_loss_n300_p100_twopercent = entropy_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_twopercent = bounded_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_twopercent = Frobenius_norm_precision_n300_p100_twopercent_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_twopercent];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_twopercent];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_twopercent];

ssIters = ssIters2 +1;
ssIters2 = ssIters ;

 end
 

 clear ssIters ssIters2 

 ssIters = 26;
 ssIters2 = ssIters + 4;

 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_twopercent_SS_%dto%d.mat', ssIters, ssIters2));

 entropy_loss_n300_p100_twopercent_tmp = entropy_loss_n300_p100_twopercent(:);
bounded_loss_n300_p100_twopercent_tmp = bounded_loss_n300_p100_twopercent(:);
Frobenius_norm_precision_n300_p100_twopercent_tmp = Frobenius_norm_precision_n300_p100_twopercent(:);

entropy_loss_n300_p100_twopercent = entropy_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_twopercent = bounded_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_twopercent = Frobenius_norm_precision_n300_p100_twopercent_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_twopercent];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_twopercent];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_twopercent];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 4;

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

table_p100_n300_twopercent_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p100_n300_twopercent_SS];


clearvars -except combine_tables

%twopercent p100 n300 VB 


 ssIters = 1;
 ssIters2 = ssIters + 24;
 
  ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];

 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_twopercent_VB_%dto%d.mat', ssIters, ssIters2));
 entropy_loss_n300_p100_twopercent_tmp = entropy_loss_n300_p100_twopercent(:);
bounded_loss_n300_p100_twopercent_tmp = bounded_loss_n300_p100_twopercent(:);
Frobenius_norm_precision_n300_p100_twopercent_tmp = Frobenius_norm_precision_n300_p100_twopercent(:);

entropy_loss_n300_p100_twopercent = entropy_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_twopercent = bounded_loss_n300_p100_twopercent_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_twopercent = Frobenius_norm_precision_n300_p100_twopercent_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_twopercent];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_twopercent];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_twopercent];

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
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_twopercent_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_twopercent_VB];


clearvars -except combine_tables

%twopercent p100 n300 HS 

load('CholeskyDecomp_p100_n300_twopercent_HS_final.mat');

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

table_p100_n300_twopercent_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_twopercent_HS];


clearvars -except combine_tables

%circle p100 n300 SS 

 ssIters = 1;
 ssIters2 = ssIters ;
 
 ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 25

load(sprintf('CholeskyDecomp_p100_n300_circle_SS_%dto%d.mat', ssIters, ssIters2));


entropy_loss_n300_p100_circle_tmp = entropy_loss_n300_p100_circle(:);
bounded_loss_n300_p100_circle_tmp = bounded_loss_n300_p100_circle(:);
Frobenius_norm_precision_n300_p100_circle_tmp = Frobenius_norm_precision_n300_p100_circle(:);

entropy_loss_n300_p100_circle = entropy_loss_n300_p100_circle_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_circle = bounded_loss_n300_p100_circle_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_circle = Frobenius_norm_precision_n300_p100_circle_tmp(ssIters:ssIters2);


ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_circle];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_circle];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_circle];

ssIters = ssIters2 +1;
ssIters2 = ssIters ;

 end
 

 clear ssIters ssIters2 

 ssIters = 26;
 ssIters2 = ssIters + 4;

 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_circle_SS_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n300_p100_circle_tmp = entropy_loss_n300_p100_circle(:);
bounded_loss_n300_p100_circle_tmp = bounded_loss_n300_p100_circle(:);
Frobenius_norm_precision_n300_p100_circle_tmp = Frobenius_norm_precision_n300_p100_circle(:);

entropy_loss_n300_p100_circle = entropy_loss_n300_p100_circle_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_circle = bounded_loss_n300_p100_circle_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_circle = Frobenius_norm_precision_n300_p100_circle_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_circle];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_circle];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_circle];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 4;

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

table_p100_n300_circle_SS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p100_n300_circle_SS];


clearvars -except combine_tables

%circle p100 n300 VB 


 ssIters = 1;
 ssIters2 = ssIters + 24;

  ELoss_value_temp = [];
 BLoss_value_temp = [];
 FrobLoss_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_circle_VB_%dto%d.mat', ssIters, ssIters2));

entropy_loss_n300_p100_circle_tmp = entropy_loss_n300_p100_circle(:);
bounded_loss_n300_p100_circle_tmp = bounded_loss_n300_p100_circle(:);
Frobenius_norm_precision_n300_p100_circle_tmp = Frobenius_norm_precision_n300_p100_circle(:);

entropy_loss_n300_p100_circle = entropy_loss_n300_p100_circle_tmp(ssIters:ssIters2);
bounded_loss_n300_p100_circle = bounded_loss_n300_p100_circle_tmp(ssIters:ssIters2);
Frobenius_norm_precision_n300_p100_circle = Frobenius_norm_precision_n300_p100_circle_tmp(ssIters:ssIters2);

ELoss_value_temp = [ELoss_value_temp ; entropy_loss_n300_p100_circle];
BLoss_value_temp = [BLoss_value_temp ; bounded_loss_n300_p100_circle];
FrobLoss_value_temp = [FrobLoss_value_temp; Frobenius_norm_precision_n300_p100_circle];

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
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_circle_VB = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_circle_VB];


clearvars -except combine_tables

%circle p100 n300 HS 

load('CholeskyDecomp_p100_n300_circle_HS_final.mat');

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

table_p100_n300_circle_HS = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_circle_HS];


clearvars -except combine_tables

%Try writing table as a csv file to read into R

writetable(combine_tables,'CholeskyDecomp_Boxplot_Loss_p100_n300.csv') 

