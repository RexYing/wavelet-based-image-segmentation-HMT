%% settings

N = 3;      % number of states
T = 8192;   % length of data sequence
L = 3;      % number of value in descrete distribution

%% random parameters

a_real  = abs(randn(N,N));      % transition matrix
a_real  = a_real./kron(ones(1,N),sum(a_real')');
% assert sum(a')=ones

b_real  = abs(randn(N,L));      % each row is discrete probability of state n
b_real  = diag(1./sum(b_real'))*b_real; % normlize

b_real  =[0.90 0.05 0.05
          0.05 0.90 0.05
          0.05 0.05 0.90];

pi_real = abs(rand(N,1));       % initial distribution
pi_real = pi_real/sum(pi_real); % normlize

%% generate random state sequence

q_real=zeros(1,T);
dp=pi_real;
for c=1:T
    q_real(c)=randdist(1,1,[1:N],dp);
    dp=a_real(q_real(c),:);
end 

%% generate observation sequence

o=zeros(1,T);
for t=1:T
    q_t=q_real(t);
    b_t=b_real(q_t,:);
    o(t)=randdist(1,1,[1:L],b_t);
end

figure(1); clf; hold on; 
hist(o,100);

%% HMM parameter estimation

fprintf(1,'forward iteration\n');
[Lalpha Lpf]=hmm_forword (a, b, pi, o);

fprintf(1,'backward iteration\n');
[Lbeta  Lpb]=hmm_backword(a, b, pi, o);

fprintf(1,'log probility calculated by  forward algorithm is %f\n',Lpf);
fprintf(1,'log probility calculated by backward algorithm is %f\n',Lpb);
fprintf(1,'(the result from forward algorithm should equal to that of backword algorithm\n');
