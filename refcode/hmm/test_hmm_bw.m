%% settings

T = 1000;   % length of data sequence

Nit=8192;    % maximal iteration
th=0.01;    % thread of iteration
Nrep=50;    % number of repetion for global maximal

N=2;
No=3;
a_real =[0.9 0.1
         0.2 0.8];
b_real =[0.4 0.4 0.2
         0.2 0.4 0.4];
pi_real=[0.3 0.7];

L=size(b_real,2);

%% generate random state sequence

fprintf(1,'generating random state sequence\n');
q_real=zeros(1,T);
dp=pi_real;
for c=1:T
    q_real(c)=randdist(1,1,[1:N],dp);
    dp=a_real(q_real(c),:);
end 

%% generate observation sequence

fprintf(1,'generating observation sequence\n');
o=zeros(1,T);
for t=1:T
    q_t=q_real(t);
    b_t=b_real(q_t,:);
    o(t)=randdist(1,1,[1:L],b_t);
end

figure(1); clf; hold on; 
hist(o,100);

%% HMM parameter estimation

% iteration
Lp_final=-inf;
for crep=1:Nrep
    crep
    
    % initial guess
    a  = abs(randn(N,N));   % transition matrix
    a  = a./kron(ones(1,N),sum(a')');

    b  = abs(randn(N,L));   % each row is discrete probability of state n
    b  = diag(1./sum(b'))*b; % normlize

    pi = abs(rand(N,1));    % initial distribution
    pi = pi/sum(pi);        % normlize

    % iteration
    Lpi= log(pi);
    Lp_old=-inf;
    for cit=1:Nit
    
        [a b pi Lp]=hmm_bw(a, b, pi, o);

        if exp(Lp-Lp_old)-1<th
            break;
        else
            fprintf(1,'relative improvement of log probability: %f\n', exp(Lp-Lp_old)-1);
            Lp_old=Lp;
        end
    end
    
    if Lp>Lp_final
        Lp_final=Lp
        a_final=a;
        b_final=b;
        pi_final=pi;
    end
    
end

%% print all
a_final
b_final
pi_final
Lp_final
