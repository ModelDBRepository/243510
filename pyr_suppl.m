% pyr_suppl.m
% Pyramidal cell model case (the case in the supplementary material):

load('EULER_PyraminalModel_4.mat')
t0=0;
tf=5000;
vT=-74.27;
IT=-1.359;
[ahat, that, gEhat, gIhat] =mainQIFestimator(v,t0,tf,dt,100,[C vE vI vT IT gL vL Iapp]);

plot_actual_vs_est_cond
