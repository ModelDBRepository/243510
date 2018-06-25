% Function to estimate the excitatory and the inhibitory conductances on
% the subthreshold regime when some ionic current of quadratic type are
% active in the subthreshold regime.
%
% Input parameters: 
%       v:  vector containing the voltage of the neuron (mV).
%       dt: time step such that v has been recovered (ms).
%       TimeW: time window where the conductances are supposed to be
%              stationary.
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
%       Program returns a Solution matrix with columns that, ahat, gEhat,
%       gIhat, respectivelly; and such that
%
%            that:   time vector for the estimated conductances
%            ahat:   the quadratic coefficient
%            gEhat:  estimated excitatory conductance (nS/cm^2)
%            gIhat:  estimated inhibitory conductance (nS/cm^2)


function Sol = QIFestimatorREC(v,t0,tf,dt,TimeW,neuronParameters)

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

% Solve the MLE moving the time window.
SWd2=floor(SampleWindow/2);
init=SWd2+1;
endt=length(v)-SWd2-1;

ahat=zeros(1,Nv-SampleWindow-1);
bhat=zeros(1,Nv-SampleWindow-1);
chat=zeros(1,Nv-SampleWindow-1);
that=zeros(1,Nv-SampleWindow-1);
j=1;
for i=init:endt
    vaux=v(i-SWd2:i+SWd2);
    [aaux,baux,caux]= MLEquadratic(vaux,dt,1,'No');
    ahat(j)=aaux;
    bhat(j)=baux;
    chat(j)=caux;
    that(j)=time(i);
    j=j+1;
end


% Estimated parameters
ahat=ahat*C;
const1 = -bhat*C-2*ahat*vT;
const2 = chat*C-ahat*vT^2+IT-Iapp;
gIhat= (const2-const1*vE)/(vI-vE);
gEhat = const1-gIhat;



Sol=[transpose(that), transpose(ahat), transpose(gEhat), transpose(gIhat)];