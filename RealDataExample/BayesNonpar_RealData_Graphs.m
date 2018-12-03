%Code to create the estimated graphs for the Real Data application
%
%Author: Jami Jackson Mulgrave
clear;

addpath('C:/Users/jnjac/Documents/MATLAB/paul-kassebaum-mathworks-circularGraph-3a7926b');

load('RealData_Frequentist_Wille.mat');
load('RealData_BDGraph_Wille.mat');

%HS Bayesian nonparanormal
load('RealData_CholeskyDecomp_HS_Wille.mat');

edge_matrix_finalanalysis_matrix_HS = double(edge_matrix_finalanalysis);

myColorMap = repmat([0 0 0], [p,1]);

%call the circular graph

figure;


H = circularGraph(edge_matrix_finalanalysis_matrix_HS, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes_HS.fig')

figs = openfig('RealData_circularGraph_Bayes_HS.fig');
   saveas(figs, 'RealData_circularGraph_Bayes_HS.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere

%Read in the SS data

load('RealData_CholeskyDecomp_SS_Wille.mat');


edge_matrix_estimate_matrix_SS = double(edge_matrix_estimate);

%call the circular graph

figure;


H = circularGraph(edge_matrix_estimate_matrix_SS, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes_SS.fig')

figs = openfig('RealData_circularGraph_Bayes_SS.fig');
   saveas(figs, 'RealData_circularGraph_Bayes_SS.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere

   %Read in the VB method
   
   load('RealData_CholeskyDecomp_VB_Wille.mat');


edge_matrix_estimate_matrix_VB = double(edge_matrix_estimate);

%call the circular graph

figure;


H = circularGraph(edge_matrix_estimate_matrix_VB, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes_VB.fig')

figs = openfig('RealData_circularGraph_Bayes_VB.fig');
   saveas(figs, 'RealData_circularGraph_Bayes_VB.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere


%Read in the frequentist STARS nonpararnormal model

figure;

stars_graph = circularGraph(edge_matrix_est_stars, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_stars.fig')

figs = openfig('RealData_circularGraph_Frequentist_stars.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_stars.png');
   

%Read in the frequentist RIC nonpararnormal model

figure;

ric_graph = circularGraph(edge_matrix_est_ric, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_ric.fig')

figs = openfig('RealData_circularGraph_Frequentist_ric.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_ric.png');

%Read in the frequentist EBIC nonpararnormal model


figure;

ebic_graph = circularGraph(edge_matrix_est_ebic, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_ebic.fig')

figs = openfig('RealData_circularGraph_Frequentist_ebic.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_ebic.png');%How many edges do we get for each?

   %BDGraph
   
   
   figure;

BDGraph_graph = circularGraph(edge_matrix_selected, 'Colormap',myColorMap);

savefig('RealData_circularGraph_BDGraph.fig')

figs = openfig('RealData_circularGraph_BDGraph.fig');
   saveas(figs, 'RealData_circularGraph_BDGraph.png');
   
   %How many different edges per graph?
   

indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 
  sum_edges_bayes_HS =  sum(edge_matrix_finalanalysis_matrix_HS(upperind) == 1);

    sum_edges_bayes_SS =  sum(edge_matrix_estimate_matrix_SS(upperind) == 1);

      sum_edges_bayes_VB =  sum(edge_matrix_estimate_matrix_VB(upperind) == 1);

    sum_edges_stars =  sum(edge_matrix_est_stars(upperind) == 1);

        sum_edges_ric =  sum(edge_matrix_est_ric(upperind) == 1);

          sum_edges_ebic =  sum(edge_matrix_est_ebic(upperind) == 1);
          
                              sum_edges_BDGraph =  sum(edge_matrix_selected(upperind) == 1);

