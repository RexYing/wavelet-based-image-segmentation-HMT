function [ result ] = fast_bilateral_filter( img, geomSigma, photoSigma, halfWindowSize )
% This is an alternative implementation that aims to speed up bilateral
% filtering.
%
% A series of cosine functions are used to approximate gaussian function
% for the range part of the bilateral filter, so that convolution of linear
% shift-invariant system may be used to speed up the process
%
% The rationale behind is that for cosine or equivalently, exponential,
% phi(f(y) - f(x)) can be expressed in terms of products of phi(f(y)) and phi(x),
% extract phi(x) out of the integral/sum and the remaining terms are suitable
% for convolution. Here phi(s) is the range function that measures
% photometric distances, and is usually gaussian.
%
% Approximation of Gaussian using cosine functions, which are also even 
% functions, is feasible, and more cosine terms yields results that are
% indistinguishable from gaussian due to eventual convergence of the series, 
% The range of image is assumed to be [0, 255], and all ranage values are 
% integers.
% 
% Paper in which the idea is proposed:
% Chaudhury, Kunal Narayan, Daniel Sage, and Michael Unser. 
% "Fast O (1) bilateral filtering using trigonometric range kernels." 
% arXiv preprint arXiv:1105.4204 (2011).
% 
%
% Rex
%

DEFAULT_HALF_WSIZE = 3;
if (nargin == 3)
    halfWindowSize = DEFAULT_HALF_WSIZE;
end

% Dynamic range
T = 256;
% Preset values of N: the raised power of cosine.
% (Raised power of cosine is essentially just $\sum c_k cos(kx)$ for some
% coefficients $c_k$ by trigo identity.
% The smaller the standard deviation of range kernel, the higher N has to
% be.
defaultN = 5;
maxN = 40;

% each component of phi(s) is cos(gamma f(s)) where s denotes range
% differences
% The value is set so that gamma*f(s) is within [-pi/2, pi/2], in which
% cosine functions are positive as does gaussian.
gamma = pi / (2 * T);
r = gamma * photoSigma; % just an experimental factor

% determine the raised power of cosine
if photoSigma > gamma^(-2)
    N = defaultN;
else
    N = min((gamma * photoSigma)^(-2), maxN);
end
fprintf('N (power of cosine raised) is chosen to be %d.\n', N);

% init
[numRows, numCols] = size(img);
result = zeros(size(img));
h = zeros(N + 1, numRows, numCols);
g = zeros(size(h));
d = zeros(size(h));

% Auxillary images for convolution h and g (they are functions of the intensities of
% original image on each pixel
wbar = waitbar(0, 'Bilateral Filter');
for n = 0: N
    waitbar(n / (3 * N), wbar);
    h(n + 1, :, :) = exp(1i * gamma * (2 * n - N) * double(img(:, :)) ...
        / (r * sqrt(N)) );
    
    g(n + 1, :, :) = squeeze(h(n + 1, :, :)) .* double(img(:, :));
    
    % coefficients
    d(n + 1, :, :) = 1 ./ h(n + 1, :, :) .* nchoosek(N, n) ./ (2 ^ N);
end

% filter with gaussian using SPATIAL information on auxillary images
% with convolution.
wSize = halfWindowSize * 2 + 1;
gaussianFilter = fspecial('gaussian', [wSize, wSize], geomSigma);
for n = 0: N
    waitbar(1/3 + n / (3 * N), wbar);
    h(n + 1, :, :) = imfilter(squeeze(h(n + 1, :, :)), gaussianFilter, ...
        'symmetric', 'same', 'conv');
    g(n + 1, :, :) = imfilter(squeeze(g(n + 1, :, :)), gaussianFilter, ...
        'symmetric', 'same', 'conv');
end

temp = zeros(size(img));
for n = 0: N
    waitbar(2/3 + (n + 1) / (3 * N), wbar);
    result = result + squeeze(d(n + 1, :, :) .* g(n + 1, :, :));
    temp = temp + squeeze(d(n + 1, :, :) .* h(n + 1, :, :));
end

result = result ./ temp;
% discard the complex part of elements in result, which theoretically
% should be 0, but numerically might not be so.
result = real(result);

close(wbar);

end

