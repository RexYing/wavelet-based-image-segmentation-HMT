function [ denoisedImg ] = denoise_img_ideal( img, convMat)
%DENOISE_IMG Denoise an image using convolution
%   img: matrix that represents pixels of an image
%   convMat: convolution matrix. Default is 3*3 with weights 1/9
%
%  Ideal filter in spatial domain
%
% Author: Rex

if nargin < 2
    convMat = ones(3, 3) / 9;
end

[numRowsImg, numColsImg] = size(img);
[numRows, numCols] = size(convMat);

% Mirror the matrix
img = [img(numRows: -1: 2, :); img];
img = [img; img(numRowsImg - 1: -1: numRowsImg - numRows + 1, :)];
img = [img(:, numCols: -1: 2), img];
img = [img, img(:, numColsImg - 1: -1: numColsImg - numCols + 1)];

% Convolution
img = conv2(img, convMat, 'same');

% De-mirror
denoisedImg = img(numRows: numRowsImg + numRows - 1, numCols: numColsImg + ...
    numCols - 1);
end

