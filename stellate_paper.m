% stellate_paper.m
% Stellate cell model case (the case in the paper):

load('EULER_StellarlModel_ORIGINAL_5_3g_2_Iapp-289.mat');
t0=0;
tf=5000;
vT=-58.7379;
IT=-9.496;
[ahat, that, gEhat, gIhat] =mainQIFestimator(v,t0,tf,dt,100,[C vE vI vT IT gL vL Iapp]);

plot_actual_vs_est_cond
