% Function to estimate the excitatory and the inhibitory conductances on
% the subthreshold regime when some ionic current of quadratic type are
% active in the subthreshold regime.
%
% Input parameters: 
%       v:  vector containing the voltage of the neuron (mV).
%       dt: time step such that v has been recovered (ms).
%       TimeW: time window where the conductances are supposed to be
%              stationary.
%       neuronParameters: vector containing [C vE vI vT IT Iapp] such that
%                         C is the membrane capacitance
%                         vE is the excitatory reversal potential
%                         vI is the inhibitory reversal potential
%                         (IT,vT) is the bifurcation point of the v-I curve
%                                 or the first non-zero value of tge f-I
%                                 curve
%                         Iapp is the applied current
%
%
% Output parameters:
%       Program returns a Solution matrix with columns gEhat and
%       gIhat, respectivelly; and such that
%
%            that:   time vector for the estimated conductances
%            gEhat:  estimated excitatory conductance (nS/cm^2)
%            gIhat:  estimated inhibitory conductance (nS/cm^2)


function Sol = QIFestimateGs(v,t0,tf,dt,TimeW,neuronParameters,aMLE)

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

% Solve the MLE moving the time window.
SWd2=floor(SampleWindow/2);
init=SWd2+1;
endt=length(v)-SWd2-1;

ahat=zeros(1,Nv-SampleWindow-1);
bhat=zeros(1,Nv-SampleWindow-1);
chat=zeros(1,Nv-SampleWindow-1);
that=zeros(1,Nv-SampleWindow-1);
j=1;
% for i=init:SampleWindow:endt
%     vaux=v(i-SWd2:i+SWd2);
%     [aaux,baux,caux]= MLEquadratic(vaux,dt,aMLE,'Yes');
%     ahat(1,j:j+SampleWindow-1)=aaux;
%     bhat(1,j:j+SampleWindow-1)=baux;
%     chat(1,j:j+SampleWindow-1)=caux;
%     that(1,j:j+SampleWindow-1)=time(i:i+SampleWindow-1);
%     j=j+SampleWindow;
% end
for i=init:endt
    vaux=v(i-SWd2:i+SWd2);
    [aaux,baux,caux]= MLEquadratic(vaux,dt,aMLE,'Yes');
    ahat(j)=aaux;
    bhat(j)=baux;
    chat(j)=caux;
    that(j)=time(i);
    j=j+1;
end

% Estimated parameters
const1 = -bhat*C-2*ahat*vT;
const2 = chat*C-ahat*vT^2+IT-Iapp;
gIhat= (const2-const1*vE)/(vI-vE);
gEhat = const1-gIhat;

Sol=[transpose(that), transpose(gEhat), transpose(gIhat)];