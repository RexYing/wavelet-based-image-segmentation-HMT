function [ w ] = form_wavelet_coef_mat( img )
%FORM_WAVELET_COEF_MAT Form the wavelet coefficient matrix with size equal
%to the original image size.
%   Detailed explanation goes here

[row, col] = size(img);

wname = 'haar';
% lowpass, horizontal detail, vertical detail, diagonal detail
%[CA, CH, CV, CD] = dwt2(img, wname, 'mode', 'sym');
lvl = log2(row);
[wc, s] = wavedec2(img, lvl, wname);

w = zeros(row, col);
for k = 1: lvl
    j = 2^(k-1); %j2 = j*j;

    % HH subband
    si=j+1; ei=2*j;
    sj=j+1; ej=2*j;
    w(si: ei, sj: ej) = detcoef2('h', wc, s, lvl - k + 1);

    % LH subband
    si=1; ei=j;
    sj=j+1; ej=2*j;
    w(si: ei, sj: ej) = detcoef2('v', wc, s, lvl - k + 1);

    % HL subband
    si=j+1; ei=2*j;
    sj=1; ej=j;
    w(si: ei, sj: ej) = detcoef2('d', wc, s, lvl - k + 1);
end

end

