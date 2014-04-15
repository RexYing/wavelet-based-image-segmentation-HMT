load('lenna_noisy');
degradedImg = x;
figure(1);
imshow(degradedImg);

figure(3);
result = imfilter(degradedImg, fspecial('gaussian', 5, 1), 'symmetric');
imshow(result);
