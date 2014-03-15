%% settings

N = 4;      % number of states
T = 8192;   % length of data sequence
No= 8;      % number of discrete outputs

%% random parameters

% transition matrix, a(i,j):=p(s_t+1=j|s_t=i)
a  = abs(randn(N,N));   
a  = diag(1./sum(a'))*a;

% initial distribution, pi(i):=p(s_1=i)
pi = abs(rand(N,1));    
pi = pi/sum(pi);

%  state to output probability, po(i,j):=p(o=j|s=i)
po = abs(randn(N,No));  
po = diag(1./sum(po'))*po;

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
    po_t=po(q_t,:);
    o(t)=randdist(1,1,[1:No],po_t);
end

%% test

[Lp q_hat]=hmm_vitdecd(a,pi,po,o);

%% print result

fprintf('optimal state sequence:\n')
for c=1:length(q_hat)
    fprintf('%d ',q_hat(c));
end
fprintf('\n')
fprintf(1,'error of estimate state sequence: %f\n', mean(q_hat~=q))
fprintf(1,'log probability of estimated state sequence: %f\n', Lp)
