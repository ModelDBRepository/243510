%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%      MLE for the QIF - with noise       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given v and t, this program estimate the parameters a,b,c on the SDE
%
%       dv = (a*v^2 + b*v + c)dt + sigma*dWt
%
% using the maximum likelihood estimator (MLE)
%
% v : voltage vector of length n
% dt : time vector of length n
% 

function [a,b,c] = MLEquadratic(v,dt,areal,aknown)
n=length(v);
v0=v(1:(n-1));
v1=v(2:n);
x4=sum(v0.^4)*dt;
x3=sum(v0.^3)*dt;
x2=sum(v0.^2)*dt;
x1=sum(v0)*dt;
x0=dt*(n-1);

xx2=sum((v1-v0).*(v0.^2));
xx1=sum((v1-v0).*v0);
xx0=sum((v1-v0));



if strcmp(aknown,'Yes')
    const1=xx1-areal*x3;
    const2=xx0-areal*x2;
    a=areal; 
    c=(x2*const2-x1*const1)/(x0*x2-x1*x1);
    b=(const1-x1*c)/x2; 
    
    
%     Asystem = [x2 x1; x1 x0];
%     bsystem = [xx1 - x3*areal;xx0 - x2*areal];   
%     Sol = linsolve(Asystem,bsystem);
%         
%     a=areal; b=Sol(1); c=Sol(2);

    
    
    
else
    % Solve the system to obtain \hat(a), \hat(b) and \hat(c)
    Asystem = [x4 x3 x2; x3 x2 x1; x2 x1 x0];
    bsystem = [xx2;xx1;xx0];
    %     Sol = linsolve(Asystem,bsystem);
    
    dA=det(Asystem);
    % dA= x4*x2*x0-x4*x1^2+2*x3*x1*x2-x3^2*x0-x2^3;     %det(Asystem);
    
    Ainv=[x0*x2-x1^2,-(x0*x3-x1*x2),x1*x3-x2^2;
        -(x0*x3-x1*x2),x0*x4-x2^2,-(x1*x4-x3*x2);
        x1*x3-x2^2,-(x1*x4-x3*x2),x2*x4-x3^2]/dA;
    
    Sol = Ainv*bsystem;
    
    a=Sol(1); b=Sol(2); c=Sol(3);
end

