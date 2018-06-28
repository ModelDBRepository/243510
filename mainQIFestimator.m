% Main program to estimate the excitatory and the inhibitory conductances on
% the subthreshold regime when some ionic current of quadratic type are
% active in the subthreshold regime.
%
% Input parameters: 
%       v:  vector containing the voltage of the neuron (mV).
%       dt: the time step that v has been recovered (ms).
%       TimeW: the time window where the conductances are supposed to be
%              stationary (100ms in our experimental data).
%       neuronParameters: vector containing [C vE vI vT IT gL vL Iapp] such that
%                         C is the membrane capacitance (constant)
%                         vE is the excitatory reversal potential (constant)
%                         vI is the inhibitory reversal potential (constant)
%                         (IT,vT) is the bifurcation point of the v-I curve
%                                 or the first non-zero value of tge f-I
%                                 curve (both are constant)
%                         gL leak conductance (constant)
%                         vL leak reversal potential (constant)
%                         Iapp is the applied current (constant)
%
%
% Output parameters:
%       Program returns a constant value ahat, and vectors that, gEhat,
%       gIhat, such that
%
%            ahat:   the quadratic coefficient
%            that:   time vector for the estimated conductances
%            gEhat:  estimated excitatory conductance (nS/cm^2)
%            gIhat:  estimated inhibitory conductance (nS/cm^2)


function [ahat, that, gEhat, gIhat] = mainQIFestimator(v,t0,tf,dt,TimeW,neuronParameters)

% Parameters of the neuron needed to estimate the conductances.
C=neuronParameters(1);
vE=neuronParameters(2);
vI=neuronParameters(3);
vT=neuronParameters(4);
IT=neuronParameters(5);
gL=neuronParameters(6);
vL=neuronParameters(7);
Iapp=neuronParameters(8);


% Time parameters and Sample Window for the MLE
Nv=length(v);
time=t0:dt:tf;
SampleWindow=TimeW/dt;


if dt <= 0
    error('The Time step must be a non-zero positive integer.')
end
if dt > 1
    warning('The Time step is to big to work with precision.')
end
if TimeW < 0
    error('The Time Window must be a positive integer.')
end
if Nv-SampleWindow-1 < 1
    error('The Time Window is too big according to the voltage data.')
end

% 1. initialize a as the mean value when the 3 parameters are estimated
% (theta=(a,b,c)):

Sol= QIFestimatorREC(v,t0,tf,dt,TimeW,neuronParameters);

that0=Sol(:,1);
ahat=abs(mean(Sol(:,2)));
gEhat0=filtrat(TimeW/dt,Sol(:,3),that0);
gIhat0=filtrat(TimeW/dt,Sol(:,4),that0);


% 2. Estimate g_E, g_I for an specific a (theta=(b,c))
Sol = QIFestimateGs(v,t0,tf,dt,TimeW,neuronParameters,ahat);
that=Sol(:,1);
gEhat=Sol(:,2);
gIhat=Sol(:,3);

% 3. Filter the results to smooth the traces...
disp('Would you like to filter your results (using a median filter)? (surround answer in single quotes)')
mfilter = input('Yes (Y), No (N): ');
if strcmp(mfilter,'Y') || strcmp(mfilter,'y') || strcmp(mfilter,'Yes') || strcmp(mfilter,'yes')
    nfilt=(TimeW/dt)/2;
    [gEhat,that]=filtrat(nfilt,gEhat,that);
    [gIhat,that]=filtrat(nfilt,gIhat,that);
end

figure();
hold on;
subplot(2,1,1)
plot(that,gEhat,'-','Color',[0.4 0.4 1],'LineWidth',2);
xlabel('time (ms)','FontSize',16);
ylabel(' g_E(t) (mS/cm^2)','FontSize',16);
set(gca,'FontSize',14);
hold off;

subplot(2,1,2)
plot(that,gIhat,'-','Color',[1,0.4,0.6],'LineWidth',2);
xlabel('time (ms)','FontSize',16);
ylabel('g_I(t) (mS/cm^2)','FontSize',16);
set(gca,'FontSize',14);
hold off;
