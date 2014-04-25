%% Create test images for segmentation

% diagonal
img1 = double(rgb2gray(imread('data/glass/glass1_t0.png')) ) / 255;
img2 = double(rgb2gray(imread('data/bark/bark13_t0.png')) ) / 255;

% assume both images are squares and 
row = size(img1, 1);
img = zeros(size(img1));
for i = 1: row
    for j = 1: row
        if (i + j < row / 2) || (i + j > 3 * row / 2)
            img(i, j) = img1(i, j);
        else
            img(i, j) = img2(i, j);
            %img(i, j) = 0.5;
        end
    end
end
imshow(img);
imwrite(img, 'data/mytest/test1.png');
