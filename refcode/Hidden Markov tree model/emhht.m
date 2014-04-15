function [ESN,PSN,MUN,SIN]=emhht(w,ES,PS,MU,SI,zm)
%function [ESN,PSN,MUN,SIN]=emhht(w,ES,PS,MU,SI,zm)
% updates HH subband of HMT model once
%
% the data structure of 2-D DWT follows the format used
% in the Rice wavelelet matlab toolbox (see README)
%
% internal variables:
% M : no. of mixture densities
% P : size of image (PxP pixels)
% level : no. of levels of HMT
%
% Input :
% w : data PxP matrix (wavelet transform of an image)
% ES : state transition matrix (MxMxPxP)
% PS : state probability matrix (MxPxP)
% MU : mean matrix (MxPxP)
% SI : variance matrix (MxPxP)
% zm : type of density functions
%  zm = 1 : zero mean (do not update MU)
%  zm = 0 : nonzero mean (update MU)
%
% Output:
% Updated ES PS MU SI in ESN PSN MUN SIN

M=size(ES,1);
P=size(w,1);
level=log2(P);
BE=zeros(M,P,P);
BEP=zeros(M,P,P);
BER=zeros(M,P,P);

AL=zeros(M,P,P);
P1=zeros(M,P,P);
P2=zeros(M,M,P,P);

%UP step
wtmp = repmat(w,[1 1 M]);
wtmp = shiftdim(wtmp,2);
si=2^(level-1)+1; ei=P; sj=2^(level-1)+1; ej=P;
gtmp = gauss(wtmp,MU,SI);
scale = repmat(mean(gtmp,1),[M 1 1]);
BE(:,si:ei,sj:ej) = gtmp(:,si:ei,sj:ej)./scale(:,si:ei,sj:ej);

%clear MUtmp SItmp;

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
    scale(:,sni:eni,snj:enj) = scale(:,sni:eni,snj:enj).*scaletmp;
    % Mixed Gaussian of the parent, times all its child beta probability
    BE(:,sni:eni,snj:enj)=gtmp(:,sni:eni,snj:enj)./scale(:,sni:eni,snj:enj).*BCtmp;
    
    %construct BE(:,pai(i),paj(j),dindex) matrix
    Btmp=zeros(M,J,J);
    % BEP (the probability of the tree rooted on one child given the parent)
    % has 4 times entries than BE in quad-tree structure. So expand BE.
    Btmp(:,1:2:J,1:2:J)=BE(:,sni:eni,snj:enj);
    Btmp(:,2:2:J,1:2:J)=BE(:,sni:eni,snj:enj);
    Btmp(:,1:2:J,2:2:J)=BE(:,sni:eni,snj:enj);
    Btmp(:,2:2:J,2:2:J)=BE(:,sni:eni,snj:enj);
    % probability of the set of coefficients obtained by removing one subtree
    % T_{ij} from its parent
    BER(:,si:ei,sj:ej)=Btmp./BEP(:,si:ei,sj:ej);
end;

clear EStmp BEtmp BCtmp Btmp;

%DOWN step
%initialize at the coarsest scale (the dc value would be PS(:, 1, 1) )
% The coarsest HH component is PS(:, 2, 2).
AL(:,2,2) = PS(:,2,2);

for k=2:level
    J = 2^(k-1); J2=J*J;
    si=J+1; ei=2*J; sj=J+1; ej=2*J;
    sni = J/2+1; eni = si-1; snj = J/2+1; enj = sj-1;
    Atmp=zeros(M,J,J);
    % duplicate entries from the next coarser scale
    Atmp(:,1:2:J,1:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,2:2:J,1:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,1:2:J,2:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,2:2:J,2:2:J)=AL(:,sni:eni,snj:enj);
    
    % multiply by BER for each state (M of them)
    Atmp = repmat(reshape(Atmp.*BER(:,si:ei,sj:ej),1,M*J2),[M 1]);
    EStmp = reshape(ES(:,:,si:ei,sj:ej),M,M*J2);
    ALtmp = reshape(EStmp.*Atmp,[M M J J]);
    AL(:,si:ei,sj:ej) = squeeze(sum(ALtmp,2));
end;

clear Atmp EStmp ALtmp;

%compute probabilities
for k=2:level
    J=2^(k-1); J2=J*J;
    si=J+1; ei=2*J; sj=J+1; ej=2*J;
    sni = J/2+1; eni = si-1; snj = J/2+1; enj = sj-1;
    
    % given the state of this node, the subtree rooted at the node and the
    % rest of the tree are independent. So multiply them together gives the
    % probability of a node at given state and the entire tree given the
    % hidden parameters.
    temp = repmat(sum(AL(:,si:ei,sj:ej).*BE(:,si:ei,sj:ej), 1),[M 1]);
    % also need to devide by the sum of all products
    % since P(\mathcal{T}_1 | \theta) is the sum of all products
    P1(:,si:ei,sj:ej) = AL(:,si:ei,sj:ej).*BE(:,si:ei,sj:ej)./temp;
    
    %compute P2
    
    % prepare for multiplication
    if M == 2
        % For M=2 the following may be faster
        BEtmp = zeros(M,M*J2);
        BEtmp(:,1:M:(M*J2))=reshape(BE(:,si:ei,sj:ej),M,J2);
        BEtmp(:,2:M:(M*J2))=BEtmp(:,1:M:(M*J2));
    else
        % For general M (not equal to 2) use the following
        BEtmp = zeros(M,M*J2);
        for m=1:M
            BEtmp(:,m:M:(M*J2))=reshape(BE(:,si:ei,sj:ej,:),M,J2);
        end;
    end;
    
    BEtmp = reshape(BEtmp,[M M J J]);
    
    EStmp = ES(:,:,si:ei,sj:ej);
    % expand alpha_{parent}
    Atmp=zeros(M,J,J);
    Atmp(:,1:2:J,1:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,2:2:J,1:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,1:2:J,2:2:J)=AL(:,sni:eni,snj:enj);
    Atmp(:,2:2:J,2:2:J)=AL(:,sni:eni,snj:enj);
    Atmp = repmat(reshape(Atmp,1,M*J2),[M 1]);
    Atmp = reshape(Atmp,[M M J J]);
    BERtmp = repmat(reshape(BER(:,si:ei,sj:ej),1,M*J2),[M 1]);
    BERtmp = reshape(BERtmp,[M M J J]);
    temp = repmat(reshape(temp,1,M*J2),[M 1]);
    temp = reshape(temp, [M M J J]);
    % the formula for 
    % P(S_i = m, S_{\rou(i)} = n, \mathcal{T}_1 | \theta)
    % similarly, need to devide by the sum of all 
    % P(S_i = m, \mathcal{T}_1 | \theta) 
    % (only 2 of them if M = 2)
    P2(:,:,si:ei,sj:ej)=BEtmp.*EStmp.*Atmp.*BERtmp./temp;
end;
% Base case
P1(:,2,2)=AL(:,2,2).*BE(:,2,2)./repmat(sum(AL(:,2,2).*BE(:,2,2),1),[M 1 1]);

clear temp BEtmp EStmp Atmp BERtmp;

%M step
PS(:,2,2)=P1(:,2,2);

% all the averaging is because of the Gaussian distribution
for k=2:level
    J=2^(k-1); J2=J*J;
    si=J+1; ei=2*J; sj=J+1; ej=2*J;
    sni = J/2+1; eni = si-1; snj = J/2+1; enj = sj-1;
    
    % update probability of states
    % at each scale sum together all possible 
    pstmp = sum(sum(P1(:,si:ei,sj:ej),3),2)/J2;
    pstmp = pstmp.*(pstmp>1e-4)+1e-4*(pstmp<=1e-4);
    PS(:,si:ei,sj:ej) = repmat(pstmp,[1 J J]);
    
    if zm == 0
        % do not update MU if zero mean densities
        mutmp = sum(sum(wtmp(:,si:ei,sj:ej).*P1(:,si:ei,sj:ej),3),2)/J2;
        MU(:,si:ei,sj:ej) = repmat(mutmp,[1 J J])./PS(:,si:ei,sj:ej);
    end;
    
    % weighting averaging of empirical variances with weight chosen (P1) in
    % proportion to the probability of each mixture.
    sitmp = sum(sum((wtmp(:,si:ei,sj:ej)-MU(:,si:ei,sj:ej)).^2.*P1(:,si:ei,sj:ej),3),2)/J2;
    % divide by the sum (which is the probability for each state)
    SI(:,si:ei,sj:ej) = repmat(sitmp,[1 J J])./PS(:,si:ei,sj:ej);
    
    estmp =sum(sum(P2(:,:,si:ei,sj:ej),4),3)/J2;
    ptmp = [PS(:,sni,snj)'; PS(:,sni,snj)'];
    ES(:,:,si:ei,sj:ej)= repmat(estmp,[1 1 J J])./repmat(ptmp,[1 1 J J]);
end; %k

ESN=ES; PSN=PS; MUN=MU; SIN=SI;








