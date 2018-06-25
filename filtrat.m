function [vfilt,tfilt]=filtrat(nfilt,v,t)

% filtrat       Median filter comput a special median filter
%
%       Inputs:
%
%               nfilt = the sample vector length
%               v = the vector which we want to pass the filter
%               t = the time vector in this case (for this reason is 
%                   special filter).
%
%       Outputs:
%
%               vfilt = the vector v filtered
%               tfilt = the corresponding time vector to vfilt


posfilt=ceil(nfilt/2);
n=length(v);
if (n>nfilt+1)
    for i=1:posfilt-1
        vfilt(i)=v(i);
        tfilt(i)=t(i);
    end
    for i=n-posfilt+1:n
        vfilt(i)=v(i);
        tfilt(i)=t(i);
    end
    for i=posfilt:n-posfilt
        tfilt(i)=t(i);
        for j=-posfilt+1:posfilt-1
            A(j+posfilt+1)=v(i+j);
        end
        [B,IX]=sort(A);
        if (posfilt ~= nfilt/2) %odd case
            vfilt(i)=B(posfilt);
        else %even case
            vfilt(i)=(B(posfilt)+B(posfilt+1))/2;
        end
    end
else
    fprintf(1,'you have to take a bigger step');
end
