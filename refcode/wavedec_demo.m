% 2D image analysis.

% Load a test image.  Matlab test images consist of a matrix, X,
% color palette, map, which maps each value of the matrix to a
% color.  Here, we will apply the Discrete Wavelet Transform to X.
load woman2
% load lena
X = x * 255;
%load detfingr; X = X(1:200,51:250);

close all
clf
image(X)
%colormap(map)
colormap('gray');
axis image; set(gca,'XTick',[],'YTick',[]); title('Original')
pause

% We will use the 9/7 filters with symmetric extension at the
% boundaries.
dwtmode('sym')
wname = 'bior4.4'

% Plot the structure of a two stage filter bank.
t = wtree(X,2,'bior4.4');
plot(t)
pause
close(2)

% Compute a 2-level decomposition of the image using the 9/7 filters.
[wc,s] = wavedec2(X,2,wname);

% Extract the level 1 coefficients.
a1 = appcoef2(wc,s,wname,1);         
h1 = detcoef2('h',wc,s,1);           
v1 = detcoef2('v',wc,s,1);           
d1 = detcoef2('d',wc,s,1);           

% Extract the level 2 coefficients.
a2 = appcoef2(wc,s,wname,2);
h2 = detcoef2('h',wc,s,2);
v2 = detcoef2('v',wc,s,2);
d2 = detcoef2('d',wc,s,2);

% Display the decomposition up to level 1 only.
ncolors = size(map,1);              % Number of colors.
sz = size(X);
cod_a1 = wcodemat(a1,ncolors); cod_a1 = wkeep(cod_a1, sz/2);
cod_h1 = wcodemat(h1,ncolors); cod_h1 = wkeep(cod_h1, sz/2);
cod_v1 = wcodemat(v1,ncolors); cod_v1 = wkeep(cod_v1, sz/2);
cod_d1 = wcodemat(d1,ncolors); cod_d1 = wkeep(cod_d1, sz/2);
image([cod_a1,cod_h1;cod_v1,cod_d1]);
axis image; set(gca,'XTick',[],'YTick',[]); title('Single stage decomposition')
colormap(map)
pause

% Display the entire decomposition upto level 2.
cod_a2 = wcodemat(a2,ncolors); cod_a2 = wkeep(cod_a2, sz/4);
cod_h2 = wcodemat(h2,ncolors); cod_h2 = wkeep(cod_h2, sz/4);
cod_v2 = wcodemat(v2,ncolors); cod_v2 = wkeep(cod_v2, sz/4);
cod_d2 = wcodemat(d2,ncolors); cod_d2 = wkeep(cod_d2, sz/4);
image([[cod_a2,cod_h2;cod_v2,cod_d2],cod_h1;cod_v1,cod_d1]);
axis image; set(gca,'XTick',[],'YTick',[]); title('Two stage decomposition')
colormap(map)
pause

% Here are the reconstructed branches
ra2 = wrcoef2('a',wc,s,wname,2);
rh2 = wrcoef2('h',wc,s,wname,2);
rv2 = wrcoef2('v',wc,s,wname,2);
rd2 = wrcoef2('d',wc,s,wname,2);

ra1 = wrcoef2('a',wc,s,wname,1);
rh1 = wrcoef2('h',wc,s,wname,1);
rv1 = wrcoef2('v',wc,s,wname,1);
rd1 = wrcoef2('d',wc,s,wname,1);

cod_ra2 = wcodemat(ra2,ncolors);
cod_rh2 = wcodemat(rh2,ncolors);
cod_rv2 = wcodemat(rv2,ncolors);
cod_rd2 = wcodemat(rd2,ncolors);
cod_ra1 = wcodemat(ra1,ncolors);
cod_rh1 = wcodemat(rh1,ncolors);
cod_rv1 = wcodemat(rv1,ncolors);
cod_rd1 = wcodemat(rd1,ncolors);
subplot(3,4,1); image(X); axis image; set(gca,'XTick',[],'YTick',[]); title('Original')
subplot(3,4,5); image(cod_ra1); axis image; set(gca,'XTick',[],'YTick',[]); title('ra1')
subplot(3,4,6); image(cod_rh1); axis image; set(gca,'XTick',[],'YTick',[]); title('rh1')
subplot(3,4,7); image(cod_rv1); axis image; set(gca,'XTick',[],'YTick',[]); title('rv1')
subplot(3,4,8); image(cod_rd1); axis image; set(gca,'XTick',[],'YTick',[]); title('rd1')
subplot(3,4,9); image(cod_ra2); axis image; set(gca,'XTick',[],'YTick',[]); title('ra2')
subplot(3,4,10); image(cod_rh2); axis image; set(gca,'XTick',[],'YTick',[]); title('rh2')
subplot(3,4,11); image(cod_rv2); axis image; set(gca,'XTick',[],'YTick',[]); title('rv2')
subplot(3,4,12); image(cod_rd2); axis image; set(gca,'XTick',[],'YTick',[]); title('rd2')
pause

% Adding together the reconstructed average at level 2 and all of
% the reconstructed details gives the full reconstructed image.
Xhat = ra2 + rh2 + rv2 + rd2 + rh1 + rv1 + rd1;
sprintf('Reconstruction error (using wrcoef2) = %g', max(max(abs(X-Xhat))))

% Another way to reconstruct the image.
XXhat = waverec2(wc,s,wname);
sprintf('Reconstruction error (using waverec2) = %g', max(max(abs(X-XXhat))))

% Compression can be accomplished by applying a threshold to the
% wavelet coefficients.  wdencmp is the function that does this.
% 'h' means use hard thresholding. Last argument = 1 means do not
% threshold the approximation coefficients.
%    perfL2 = energy recovery = 100 * ||wc_comp||^2 / ||wc||^2.
%             ||.|| is the L2 vector norm.
%    perf0 = compression performance = Percentage of zeros in wc_comp.
thr = 20;                                                    
[X_comp,wc_comp,s_comp,perf0,perfL2] = wdencmp('gbl',wc,s,wname,2,thr,'h',1);

%% Denoise
[thr,sorh,keepapp] = ddencmp('den','wv',X);
xd = wdencmp('gbl',X,'sym4',2,thr,sorh,keepapp);
imshow(xd, [0, 255]);

%%
clf
subplot(1,2,1); image(X); axis image; set(gca,'XTick',[],'YTick',[]);
title('Original')
cod_X_comp = wcodemat(X_comp,ncolors);
subplot(1,2,2); image(cod_X_comp); axis image; set(gca,'XTick',[],'YTick',[]);
title('Compressed using global hard threshold')
xlabel(sprintf('Energy retained = %2.1f%% \nNull coefficients = %2.1f%%',perfL2,perf0))
pause

% Better compression can be often be obtained if different thresholds
% are allowed for different subbands.
thr_h = [21 17];        % horizontal thresholds.              
thr_d = [23 19];        % diagonal thresholds.                
thr_v = [21 17];        % vertical thresholds.                
thr = [thr_h; thr_d; thr_v];
% alternative: the default value by ddencmp
% [thr,sorh,keepapp] = ddencmp('den','wv',X);
[X_comp,wc_comp,s_comp,perf0,perfL2] = wdencmp('lvd',X,wname,2,thr,'h');

clf
subplot(1,2,1); image(X); axis image; set(gca,'XTick',[],'YTick',[]);
title('Original')
cod_X_comp = wcodemat(X_comp,ncolors);
subplot(1,2,2); image(cod_X_comp); axis image; set(gca,'XTick',[],'YTick',[]);
title('Compressed using variable hard thresholds')
xlabel(sprintf('Energy retained = %2.1f%% \nNull coefficients = %2.1f%%',perfL2,perf0))

% Return to default settings.
dwtmode('zpd')
