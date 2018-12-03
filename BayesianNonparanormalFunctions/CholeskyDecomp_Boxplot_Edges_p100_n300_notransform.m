%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the CholeskyDecomp paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AR2 p100 n300 SS 


 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 SP_value_temp = [];
 SE_value_temp = [];
 MCC_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_AR2_SS_notransform_%dto%d.mat', ssIters, ssIters2));

SP_matrix = SP_matrix(:,:,ssIters:ssIters2);
SE_matrix = SE_matrix(:,:,ssIters:ssIters2);
MCC_matrix = MCC_matrix(:,:,ssIters:ssIters2);

SP_value_temp = [SP_value_temp; SP_matrix];
SE_value_temp = [SE_value_temp; SE_matrix];
MCC_value_temp =[MCC_value_temp; MCC_matrix];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 

SP_value = reshape(SP_value_temp, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_value_temp, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_value_temp, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_AR2_SS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [table_p100_n300_AR2_SS];


clearvars -except combine_tables

%AR2 p100 n300 VB 

load('CholeskyDecomp_p100_n300_AR2_VB_notransform.mat');

SP_value = reshape(SP_matrix, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_matrix, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_matrix, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);


Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_AR2_VB = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_AR2_VB];


clearvars -except combine_tables

%AR2 p100 n300 HS 

load('CholeskyDecomp_p100_n300_AR2_HS_notransform_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR2'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_AR2_HS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_AR2_HS];


clearvars -except combine_tables

%twopercent p100 n300 SS 


 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 SP_value_temp = [];
 SE_value_temp = [];
 MCC_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_twopercent_SS_notransform_%dto%d.mat', ssIters, ssIters2));

SP_matrix = SP_matrix(:,:,ssIters:ssIters2);
SE_matrix = SE_matrix(:,:,ssIters:ssIters2);
MCC_matrix = MCC_matrix(:,:,ssIters:ssIters2);

SP_value_temp = [SP_value_temp; SP_matrix];
SE_value_temp = [SE_value_temp; SE_matrix];
MCC_value_temp =[MCC_value_temp; MCC_matrix];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 

SP_value = reshape(SP_value_temp, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_value_temp, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_value_temp, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_twopercent_SS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables ; table_p100_n300_twopercent_SS];


clearvars -except combine_tables

%twopercent p100 n300 VB 

load('CholeskyDecomp_p100_n300_twopercent_VB_notransform.mat');

SP_value = reshape(SP_matrix, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_matrix, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_matrix, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);


Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_twopercent_VB = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_twopercent_VB];


clearvars -except combine_tables

%twopercent p100 n300 HS 

load('CholeskyDecomp_p100_n300_twopercent_HS_notransform_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_twopercent_HS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_twopercent_HS];


clearvars -except combine_tables

%circle p100 n300 SS 


 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 SP_value_temp = [];
 SE_value_temp = [];
 MCC_value_temp = [];
 
 while ssIters <= 100

load(sprintf('CholeskyDecomp_p100_n300_circle_SS_notransform_%dto%d.mat', ssIters, ssIters2));

SP_matrix = SP_matrix(:,:,ssIters:ssIters2);
SE_matrix = SE_matrix(:,:,ssIters:ssIters2);
MCC_matrix = MCC_matrix(:,:,ssIters:ssIters2);

SP_value_temp = [SP_value_temp; SP_matrix];
SE_value_temp = [SE_value_temp; SE_matrix];
MCC_value_temp =[MCC_value_temp; MCC_matrix];

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 

SP_value = reshape(SP_value_temp, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_value_temp, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_value_temp, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Bernoulli-Gaussian'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_circle_SS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_circle_SS];


clearvars -except combine_tables

%circle p100 n300 VB 

load('CholeskyDecomp_p100_n300_circle_VB_notransform.mat');

SP_value = reshape(SP_matrix, [reps,1]);
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = reshape(SE_matrix, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = reshape(MCC_matrix, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);


Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Variational Bayes'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_circle_VB = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_circle_VB];


clearvars -except combine_tables

%circle p100 n300 HS 

load('CholeskyDecomp_p100_n300_circle_HS_notransform_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Horseshoe'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n300_circle_HS = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n300_circle_HS];


clearvars -except combine_tables

%Try writing table as a csv file to read into R

writetable(combine_tables,'CholeskyDecomp_Boxplot_Edges_p100_n300_notransform.csv') 

