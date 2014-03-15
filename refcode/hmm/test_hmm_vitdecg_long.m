%% settings

N = 4;      % number of states
T = 8192;  % length of data sequence

%% random parameters

a  = abs(randn(N,N));   % transition matrix
a  = a./kron(ones(1,N),sum(a')');
% assert sum(a')=ones

mu = 5*randn(N,1);        % mean of Guassian at 4 states
var= abs(randn(N,1));   % variance of Gaussian at 4 states

pi = abs(rand(N,1));    % initial distribution
pi = pi/sum(pi);

%% generate random state sequence

q=zeros(1,T);
dp=pi;
for c=1:T
    q(c)=randdist(1,1,[1:N],dp);
    dp=a(q(c),:);
end

%% generate observation sequence

o=zeros(1,T);
for t=1:T
    q_t=q(t);
    var_t=var(q_t);
    mu_t=mu(q_t);
    o(t)=randn()*var_t+mu_t;
end

figure(1); clf; hold on; 
hist(o,100);

%% test
[Lp,q_hat]=hmm_vitdecg(a,pi,mu,var,o);

%% print result

figure(2); clf; hold on;
plot(q(400:500), 'r.-');
plot(q_hat(400:500), 'k.-');

fprintf(1,'state estimation error rate by viterbi algorithm is %f\n', mean(q~=q_hat))

%% simple method to find state, no state transition information is used
q_hats=zeros(1,T);
for t=1:T
    dd=-log(sqrt(2*3.1415926)*var) - 0.5*(o(t)-mu).^2./var;
    [mx idx]=max(dd);
    q_hats(t)=idx;
end

fprintf(1,'state estimation error rate by simple hard decision is %f\n', mean(q~=q_hats))

% compare simple method with viterbi method (should return 1 below)
mean(q~=q_hat)<mean(q~=q_hats)
