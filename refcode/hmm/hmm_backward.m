function [Lbeta, Lp]=hmm_backword(a, b, pi, o)
%--------------------------------------------------------------------------
%Backword algorithm
%
%   [Lbeta,q] = hmm_backword(a,b,pi,o)
%
%   inputs:
%       a(i,j)  transition probability matrix, a(i,j) :=p(q_t=j|q_t-1=i)
%       pi(i)   initial probability,           pi(i)  :=p(q_1=i        )
%       b(i,j)  output probability matrix,     b(i,j) :=p(  o=j|q=i    )
%       o       observation sequence
%
%
%   outputs:
%       Lbeta   log probability of beta
%       Lp      log probability of observation sequence
%--------------------------------------------------------------------------

%% check inputs
% number of states
N=size(a,1);
if N~=size(a,2)
    fprintf(1,'error, state transition probability matrix should be square\n');
    return;
end

% number of type of discrete outputs  
No=size(b,2);
if size(b,1)~=N
    fprintf(1,'error, row size of po should equal to the number of state\n');
    return;
end

% length of observation
T=length(o);

pi=pi(:);   % make sure it is column vector

%% backword algorithm

Lbeta=zeros(N,T);
%Lbeta(:,T)=0;
for t=T-1:-1:1
    for i=1:N
        Lbeta(i,t)=log_sum(log(a(i,:)')+log(b(:,o(t+1)))+Lbeta(:,t+1));
    end
end

Lp=log_sum(Lbeta(:,1)+log(pi)+log(b(:,o(1))));


return;