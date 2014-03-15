
function [a, b, pi, Lp]=hmm_bw(a0, b0, pi0, o)

%--------------------------------------------------------------------------
%Baum Welch algorithm
%
%   [a, b, pi, Lp]=hmm_bw(a0, b0, pi0, o)
%
%   inputs:
%       a0(i,j) initial guess of transition probability matrix, a0(i,j) :=p(q_t=j|q_t-1=i)
%       b0(i,j) initial guess of output probability matrix,     b0(i,j) :=p(  o=j|q=i    )
%       pi0(i)  initial guess of initial probability,           pi0(i)  :=p(q_1=i        )
%       o       observation sequence
%
%   outputs:
%       a   updated transition probability matrix
%       b   updated output probability matrix
%       pi	updated initial probability
%       Lp  log probability of observation sequence
%--------------------------------------------------------------------------

%% check inputs

a=a0;
b=b0;
pi=pi0;

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

% number of types of observation variables
L=size(b,2);

pi=pi(:);   % make sure it is column vector

%% HMM parameter estimation

Lpi= log(pi);

% forward iteration
[Lalpha Lpf]=hmm_forword (a, b, pi, o);
    
% backward iteration
[Lbeta  Lpb]=hmm_backword(a, b, pi, o);

% evaluate gamma
Lgamma=zeros(N,T);
for t=1:T
    tmp=log_sum(Lalpha(:,t)+Lbeta(:,t));
    Lgamma(:,t)=Lalpha(:,t)+Lbeta(:,t)-tmp;
end

% evaluate kc
Lkc   =zeros(N,N,T-1);
for t=1:T-1
    for i=1:N
        for j=1:N
            Lkc(i,j,t)=Lgamma(i,t)+log(a(i,j))+log(b(j,o(t+1)))+Lbeta(j,t+1)-Lbeta(i,t);
        end
    end
end

% estimate pi
Lpi=Lgamma(:,1);
pi=exp(Lpi);

% estimate a
for i=1:N
    for j=1:N
        a(i,j)=sum(exp(Lkc(i,j,1:T-1)))/sum(exp(Lgamma(i,1:T-1)));
    end
end

% estimate b
for i=1:N
    for l=1:L
        b(i,l)=(o==l)*exp(Lgamma(i,:))'/sum(exp(Lgamma(i,:)));
    end
end

% obtain Lp
[Lalpha Lp]=hmm_forword (a, b, pi, o);

return;
