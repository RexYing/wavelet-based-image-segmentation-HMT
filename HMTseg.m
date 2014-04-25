%% HMT segmentation


%% Texture
ES = cell(3, 1); PS = cell(3, 1); MU = cell(3, 1); SI = cell(3, 1);

% img = double(rgb2gray(imread('data/man-made/mat2_t0.png')) ) / 255;
% w = form_wavelet_coef_mat(img);
% [ES{1},PS{1},MU{1},SI{1}] = hmttrain(w, 2);
% 
% img = double(rgb2gray(imread('data/textile/cloth29_t0.png')) ) / 255;
% w = form_wavelet_coef_mat(img);
% [ES{2},PS{2},MU{2},SI{2}] = hmttrain(w, 2);
% 
% img = double(rgb2gray(imread('data/man-made/roofTiles2_t0.png')) ) / 255;
% w = form_wavelet_coef_mat(img);
% [ES{3},PS{3},MU{3},SI{3}] = hmttrain(w, 2);

img1 = double(rgb2gray(imread('data/glass/glass1_t0.png')) ) / 255;
w = form_wavelet_coef_mat(img1);
[ES{1},PS{1},MU{1},SI{1}] = hmttrain(w, 2);

img2 = double(rgb2gray(imread('data/bark/bark13_t0.png')) ) / 255;
w = form_wavelet_coef_mat(img2);
[ES{2},PS{2},MU{2},SI{2}] = hmttrain(w, 2);

%% Compute likelihood
% number of types of training texture
num = 2;

%img = double(imread('data/tm1_1_1.png') ) / 255;
img = double(imread('data/mytest/test1.png') ) / 255;
[row, col] = size(img);
w1 = form_wavelet_coef_mat(img);
k = 4;
p = zeros(num, 2^(k-1), 2^(k-1));

P1(i) = zeros(2,row,col);
for i = 1: num
    
%     ptree = compute_likelihood( w1, ES{i}, PS{i}, MU{i}, SI{i} );
%     
%     J=2^(k-1); J2=J*J; si = J+1; ei = 2*J; sj = J+1; ej = 2*J;
%     p(i, :, :) = ptree(si: ei, sj: ej);
% 
%     J=2^(k-1); J2=J*J; si = J+1; ei = 2*J; sj = 1; ej = J;
%     p(i, :, :) = squeeze(p(i, :, :)) .* ptree(si: ei, sj: ej);
% 
%     J=2^(k-1); J2=J*J; si = 1; ei = J; sj = J+1; ej = 2*J;
%     p(i, :, :) = squeeze(p(i, :, :)) .* ptree(si: ei, sj: ej);
%     
%     squeeze(p(i, :, :))
	P1(i)=posthh(w,ES,PS,MU,SI, squeeze(P1(1, :, :)) );
    P1(i)=postlh(w,ES,PS,MU,SI, squeeze(P1(1, :, :)) );
    P1(i)=posthl(w,ES,PS,MU,SI, squeeze(P1(1, :, :)) );
    
end

[~, ind] = max(p, [], 1);
ind = squeeze(ind);
figure(2);
imagesc(ind);

%tmp =squeeze(ptree(3, :, :));

%% analyze wavelet coefficients
% h = cell([1, lvl]);
% v = cell([1, lvl]);
% d = cell([1, lvl]);
% for i = 1: lvl
%     h{i} = detcoef2('h',wc,s,i);
%     v{i} = detcoef2('v',wc,s,i);
%     d{i} = detcoef2('d',wc,s,i);
% end
% 
% hist(h{1}, 100);
%hist(v{lvl}, 100);
%hist(d{lvl}, 100);