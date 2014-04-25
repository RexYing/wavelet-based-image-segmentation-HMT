function [ ptree ] = compute_likelihood( w, ES, PS, MU, SI )
%COMPUTE_LIKELIHOOD Compute likelihood of the coefficients in the wavelet
%tree given the hidden parameters
%   Detailed explanation goes here
%

M = 2;
level = log2(size(w, 1));
[row, col] = size(w);

wtmp = repmat(w, [1 1 M]);
wtmp = shiftdim(wtmp,2);
gtmp = gauss(wtmp,MU,SI);
scale = repmat(mean(gtmp,1),[M 1 1]);

%% HH
si=2^(level-1)+1; ei=row; sj=2^(level-1)+1; ej=col;
BE(:,si:ei,sj:ej) = gtmp(:,si:ei,sj:ej)./scale(:,si:ei,sj:ej);

for k=level:-1:2
    J=2^(k-1);J2=J*J; si = J+1; ei = 2*J; sj = J+1; ej = 2*J;
    
    % reshape the left stochastic matrices for easy multiplication
    EStmp = reshape(ES(:,:,si:ei,sj:ej),M,M*J2);
    
    % prepare for multiplication
    if M == 2
        %%%%%% For M=2 the following is faster
        BEtmp = zeros(M,M*J2);
        BEtmp(:,1:M:(M*J2))=reshape(BE(:,si:ei,sj:ej),M,J2);
        BEtmp(:,2:M:(M*J2))=BEtmp(:,1:M:(M*J2));
    else
        % For general M (not equal to 2) use the following
        %
        BEtmp = zeros(M,M*4^(k-1));
        for m=1:M
            BEtmp(:,m:M:(M*4^(k-1)))=reshape(BE(:,si:ei,sj:ej,:),M,4^(k-1));
        end;
    end;
    
    BEtmp = reshape(EStmp.*BEtmp,[M M J J]);
    % probability of a tree rooted on parent, given parent is state i (i = 1, 2)
    BEP(:,si:ei,sj:ej) = squeeze(sum(BEtmp,1));
    
    sni = J/2+1; eni = si-1; snj = J/2+1; enj = sj-1;
    
    %construct betachild matrix
    BCtmp = BEP(:, si:2:ei, sj:2:ej);
    BCtmp = BCtmp.*BEP(:, si+1:2:ei, sj:2:ej);
    BCtmp = BCtmp.*BEP(:, si:2:ei, sj+1:2:ej);
    BCtmp = BCtmp.*BEP(:, si+1:2:ei, sj+1:2:ej);
    scaletmp = repmat(mean(BCtmp,1),[M 1 1]);
    scale1(:,sni:eni,snj:enj) = scale(:,sni:eni,snj:enj).*scaletmp;
    % Mixed Gaussian of the parent, times all its child beta probability
    BE(:,sni:eni,snj:enj)=gtmp(:,sni:eni,snj:enj)./scale1(:,sni:eni,snj:enj).*BCtmp;
end;

%% HL

si=2^(level-1)+1; ei=row; sj=1; ej=2^(level-1);
BE(:,si:ei,sj:ej) = gtmp(:,si:ei,sj:ej)./scale(:,si:ei,sj:ej);

for k=level:-1:2
    J=2^(k-1);J2=J*J; si = J+1; ei = 2*J; sj = 1; ej = J;
    
    EStmp = reshape(ES(:,:,si:ei,sj:ej),M,M*J2);
    
    if M == 2
        %%%%%% For M=2 the following is faster
        BEtmp = zeros(M,M*J2);
        BEtmp(:,1:M:(M*J2))=reshape(BE(:,si:ei,sj:ej),M,J2);
        BEtmp(:,2:M:(M*J2))=BEtmp(:,1:M:(M*J2));
    else
        % For general M (not equal to 2) use the following
        %
        BEtmp = zeros(M,M*4^(k-1));
        for m=1:M
            BEtmp(:,m:M:(M*4^(k-1)))=reshape(BE(:,si:ei,sj:ej,:),M,4^(k-1));
        end;
    end;
    
    BEtmp = reshape(EStmp.*BEtmp,[M M J J]);
    BEP(:,si:ei,sj:ej) = squeeze(sum(BEtmp,1));
    
    sni = J/2+1; eni = J; snj = 1; enj = J/2;
    
    %construct betachild matrix here
    BCtmp = BEP(:,si:2:ei,sj:2:ej);
    BCtmp = BCtmp.*BEP(:,si+1:2:ei,sj:2:ej);
    BCtmp = BCtmp.*BEP(:,si:2:ei,sj+1:2:ej);
    BCtmp = BCtmp.*BEP(:,si+1:2:ei,sj+1:2:ej);
    scaletmp = repmat(mean(BCtmp,1),[M 1 1]);
    scale1(:,sni:eni,snj:enj) = scale(:,sni:eni,snj:enj).*scaletmp;
    BE(:,sni:eni,snj:enj)=gtmp(:,sni:eni,snj:enj)./scale1(:,sni:eni,snj:enj).*BCtmp;
end;

%% HL

si=1; ei=2^(level-1); sj=2^(level-1)+1; ej=col;
BE(:,si:ei,sj:ej) = gtmp(:,si:ei,sj:ej)./scale(:,si:ei,sj:ej);

for k=level:-1:2
    J=2^(k-1);J2=J*J; si = 1; ei = J; sj = J+1; ej = 2*J;
    
    EStmp = reshape(ES(:,:,si:ei,sj:ej),M,M*J2);
    
    if M == 2
        %%%%%% For M=2 the following is faster
        BEtmp = zeros(M,M*J2);
        BEtmp(:,1:M:(M*J2))=reshape(BE(:,si:ei,sj:ej),M,J2);
        BEtmp(:,2:M:(M*J2))=BEtmp(:,1:M:(M*J2));
    else
        % For general M (not equal to 2) use the following
        %
        BEtmp = zeros(M,M*4^(k-1));
        for m=1:M
            BEtmp(:,m:M:(M*4^(k-1)))=reshape(BE(:,si:ei,sj:ej,:),M,4^(k-1));
        end;
    end;
    
    BEtmp = reshape(EStmp.*BEtmp,[M M J J]);
    BEP(:,si:ei,sj:ej) = squeeze(sum(BEtmp,1));
    
    sni = 1; eni = J/2; snj = J/2+1; enj = J;
    
    %construct betachild matrix here
    BCtmp = BEP(:,si:2:ei,sj:2:ej);
    BCtmp = BCtmp.*BEP(:,si+1:2:ei,sj:2:ej);
    BCtmp = BCtmp.*BEP(:,si:2:ei,sj+1:2:ej);
    BCtmp = BCtmp.*BEP(:,si+1:2:ei,sj+1:2:ej);
    scaletmp = repmat(mean(BCtmp,1),[M 1 1]);
    scale1(:,sni:eni,snj:enj) = scale(:,sni:eni,snj:enj).*scaletmp;
    BE(:,sni:eni,snj:enj)=gtmp(:,sni:eni,snj:enj)./scale1(:,sni:eni,snj:enj).*BCtmp;
end;

%% Compute likelihood
% P(T_i | \theta) = P(T_i | S_i = m, \theta) * P(S_i = m | \theta)

ptree = BE .* PS;
ptree = squeeze(sum(ptree, 1));

end

