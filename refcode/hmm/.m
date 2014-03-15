
function [Lp q_dec]=hmm_vitdecg(a, pi, mu, var, o)
%--------------------------------------------------------------------------
%Viterbi decoder (with gaussian observations)
%   [Lp,q] = hmm_vitdecd(a,pi,po,o)
%
%   inputs:
%       a(i,j)  transition probability matrix, a(i,j) :=p(q_t=j|q_t-1=i)
%       pi(i)   initial probability,           pi(i)  :=p(q_1=i        )
%       mu(i)   mean of gaussian probability at state i
%       var(i)  variance of gaussian probability at state i
%       o       observation sequence
%
%   outputs:
%       Lp      log probability of optimal path
%       q_dec   optimal state sequence
%--------------------------------------------------------------------------
%   See also HMM_DECD

%   $Revision: 1.0 $  $Date: 2008/10/04 00:00:00 $

%% check inputs

% number of states
N=size(a,1);
if N~=size(a,2)
    fprintf(1,'error, state transition probability matrix should be square\n');
    return;
end

% length of observation
T=length(o);

%% viterbi algorithm for state estimation

Ld=zeros(N,T);
f=zeros(N,T);

coff=-log(3.141592653589793)/2; % const

% init
Lb= coff-log(var)-0.5*(o(1)-mu).^2./var;
Ld(:,1)=log(pi)+Lb;
f(:,1)=0;

% interation
for t=2:T
    for j=1:N
        dd=Ld(:,t-1)+log(a(:,j));   % log probability of path from state 1~N to state j
        [mx idx]=max(dd);           % find optimal path to state j
        f(j,t)=idx;                 % record the optimal path to j @ time t
        
        o_t=o(t);
        mu_j=mu(j);
        var_j=var(j);
        Lb_j=coff-log(var_j)-0.5*(o_t-mu_j)^2/var_j; % log probability of output o_t given current state j
        Ld(j,t)=mx + Lb_j;          % log probability of optimal path to j with output o_t
    end
end

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
