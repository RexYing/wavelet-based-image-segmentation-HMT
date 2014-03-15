
function [Lp q_dec]=hmm_vitdecd(a, pi, po, o)
%--------------------------------------------------------------------------
%Viterbi decoder (with discrete observations)
%
%   [Lp,q] = hmm_vitdecd(a,pi,po,o)
%
%   inputs:
%       a(i,j)  transition probability matrix, a(i,j) :=p(q_t=j|q_t-1=i)
%       pi(i)   initial probability,           pi(i)  :=p(q_1=i        )
%       po(i,j) output probability matrix,     po(i,j):=p(  o=j|q=i    )
%       o       observation sequence
%
%
%   outputs:
%       Lp      log probability of optimal path
%       q_dec   optimal state sequence
%--------------------------------------------------------------------------
%   See also HMM_DECG

%   $Revision: 1.0 $  $Date: 2008/10/04 00:00:00 $

%% check inputs

% number of states
N=size(a,1);
if N~=size(a,2)
    fprintf(1,'error, state transition probability matrix should be square\n');
    return;
end

% number of type of discrete outputs  
No=size(po,2);
if size(po,1)~=N
    fprintf(1,'error, row size of po should equal to the number of state\n');
    return;
end

% length of observation
T=length(o);

%% viterbi algorithm for state estimation

Ld=zeros(N,T);      % store log probability
f=zeros(N,T);

% init
Lb=log(po(:,o(1)));     % Lb(i): the log probability of 1st output given 1st state i
Ld(:,1)=log(pi(:))+Lb;  % Ld(i,1): the log probability that 1st state is i and 1st output is o_1

% interation
for t=2:T
    for j=1:N
        dd=Ld(:,t-1)+log(a(:,j));   % log probability of path from state 1~N to state j
        [mx idx]=max(dd);           % find optimal path to state j
        f(j,t)=idx;                 % record the optimal path to j @ time t
        
        o_t=o(t);
        b_t=po(j,o_t);              % probability of output o_t given current state j
        
        Ld(j,t)=mx + log(b_t);      % log probability of optimal path to j with output o_t
    end
end

% final
[mx idx]=max(Ld(:,T));
Lp=mx;      % log probability of optimal path
q_last=idx; % optimal last state

%% arrange result
q_dec=zeros(1,T);
for t=T:-1:1
    q_dec(t)=q_last;
    q_last=f(q_last,t);
end

return;
