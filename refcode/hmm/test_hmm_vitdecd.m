%% settings

N=2;
No=3;
o=[2 2 3 1 1 2 2 2 1 3];
T=length(o);
a=[0.7 0.3
   0.4 0.6];
po=[0.1 0.4 0.5
    0.6 0.3 0.1];
pi=[0.6 0.4];

%% test

[Lp q_hat]=hmm_vitdecd(a,pi,po,o);

fprintf('optimal state sequence:\n')
for c=1:length(q_hat)
    fprintf('%d ',q_hat(c));
end
fprintf('\n')
% should print 1 1 1 2 2 1 1 1 2 1

fprintf(1,'log probability of estimated state sequence: %f\n', Lp)
