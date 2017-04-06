%% Q1. Matching
% close all; clear all
% 1) Manual
img = im2double(imread('images.jpg'));
% img = im2double(imread('test1.png'));
% img = im2double(imread('test2.png'));
imshow(img);
[X1,Y1] = ginput(5); % get interest points in imgA
[X2,Y2] = ginput(5); % get interest points in imgB
%% 2) Automatic
% a) Harris points detector
%
% b) calculate local feature descriptor from interest points

% Get image gradient magnitudes and directions at each pixel

%  Convert gradient directions into corresponding bin from the 8 bins
%  (0:45:360)

% For each interest point, iterate over each cell in 3x3 grid (??? change
% to 32x32). For each cell, extract subarray of gradient magnitudes and
% directions. Calculate bin count weighted by gradient magnitudes for each
% cell. Concatente and normalize to get features for interest point

%  Return descriptors for each interest point
%
% c) NN matching between descriptors from two different images

% For each interest point, get ratio of distance of the nearest neighbor to
% the point to the distance of the next nearest neighbor to the point. This
% involves computing the distance to the two nearest neighbors. We use the
% kd-tree implementation in Matlab to compute these distance for faster
% computation.

% Sort nearest neighbors in ascending order according to value of Lowe's ratio

% Return best match for each interest point along with confidences in order
% from most confident to least confident

%% Use self-written function
patch_size = 32;
imgA = 'img1.pgm';
imgB = 'img2.pgm';
threshold = 0.97;

ptA = get_interePt(imgA, patch_size);
ptB = get_interePt(imgB, patch_size);
featuresA = get_feature(imgA,ptA,patch_size);
featuresB = get_feature(imgB,ptB,patch_size);
[matchmy, confidence,dist,r] = knn_match(featuresA, featuresB, threshold);

% for match accuracy test
point = 1; % change this number for any point we want
a=ptA(:,matchmy(point,1))';
figure(1);hold on;plot(a(1,1),a(1,2),'r*')
b=ptB(:,matchmy(point,2))';
figure(2);hold on;plot(b(1),b(2),'r*')

%% Use builtin function
imgA = im2double(imread(imgA));
cornersA = detectHarrisFeatures(imgA);
figure;imshow(imgA); hold on;plot(cornersA.selectStrongest(50));

imgB = im2double(imread(imgB));
cornersB = detectHarrisFeatures(imgB);
figure;imshow(imgB); hold on;plot(cornersB.selectStrongest(50));

I1 = imgA;
I2 = imgB;
points1 = detectHarrisFeatures(I1);
points2 = detectHarrisFeatures(I2);
[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);
indexPairs = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);

%% 3) Transformation estimation

% a) homography matrix

% Assume we have two images '1' and '2' and four feature pairs:
% (p1 is the reference image)
% corre1 is a 2*n matrix. Each row is the coordinates of a point in '1'
% corre2 is a 2*n matrix. Each row is the coordinates of a point in '2'
% H is the homography matrix, such that:
% corre1_homogeneous = H * [corre2 ones(size(corre2, 1), 1)]'
n = size(corre1, 2);
if n < 3
 error('Not enough points');
end 

% the original transformation formula is:
% img2 = h*img1 where H is 3x3 matrix
% for the convenience of calculation expand H into col vector and reshape X
% after reshape and expand : x = X*H where H is 9x1 vector

X = zeros(n*3,9); % reshape of img1 denoted as X
x = zeros(n*3,1); % corresponding pt in img2 denoted as x
for i=1:n 
 X(3*(i-1)+1,1:3) = [corre2(:,i)',1];
 X(3*(i-1)+2,4:6) = [corre2(:,i)',1];
 X(3*(i-1)+3,7:9) = [corre2(:,i)',1];
 x(3*(i-1)+1:3*(i-1)+3) = [corre1(:,i);1];
end
H = (X\x)'; % obtained from x = X*H
H = reshape(H,3,3)';
% for transformation test:
i = 1; % change this number to any point you want to test
corre = H*[corre1(:,i);1];

%% by using SVD solution to compute H
% cannot get corresponding point when apply img2 = H*img1
H = get_homography(corre1,corre2);
% for transformation test:
i = 1; % change this number to any point you want to test
corre = H*[corre1(:,i);1];


