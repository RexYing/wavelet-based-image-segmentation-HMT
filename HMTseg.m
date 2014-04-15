%% Texture
img = double(rgb2gray(imread('data/textures/mat2.jpg')) );
[row, col] = size(img);
img = img(3: row - 2, 3: col - 2);
imshow(img, [0 255]);

wname = 'haar';
% lowpass, horizontal detail, vertical detail, diagonal detail
%[CA, CH, CV, CD] = dwt2(img, wname, 'mode', 'sym');
lvl = 5;
[wc,s] = wavedec2(img, lvl, wname);

%% analyze wavelet coefficients
h = cell([1, lvl]);
v = cell([1, lvl]);
d = cell([1, lvl]);
for i = 1: lvl
    h{i} = detcoef2('h',wc,s,i);
    v{i} = detcoef2('v',wc,s,i);
    d{i} = detcoef2('d',wc,s,i);
end

hist(h{1}, 100);
%hist(v{lvl}, 100);
%hist(d{lvl}, 100);