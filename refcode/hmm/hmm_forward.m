
function [Lalpha, Lp]=hmm_forword(a, b, pi, o)
%--------------------------------------------------------------------------
%Backword algorithm
%
%   [Lalpha, Lp] = hmm_forword(a,b,pi,o)
%
%   inputs:
%       a(i,j)  transition probability matrix, a(i,j) :=p(q_t=j|q_t-1=i)
%       pi(i)   initial probability,           pi(i)  :=p(q_1=i        )
%       b(i,j)  output probability matrix,     b(i,j) :=p(  o=j|q=i    )
%       o       observation sequence
%
%
%   outputs:
%       Lalpha  log probability of alpha
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

%% forword algorithm

Lalpha=zeros(N,T);
Lalpha(:,1)=log(pi)+log(b(:,o(1)));
for t=1:T-1
    for j=1:N
        Lalpha(j,t+1)=log_sum(Lalpha(:,t)+log(a(:,j)))+log(b(j,o(t+1)));
    end
end

Lp=log_sum(Lalpha(:,T));

return;
