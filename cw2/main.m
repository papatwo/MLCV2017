%% Q1. Matching
% close all; clear all
% 1) Manual

%% 2) Automatic
% a) Harris points detector

%%
% b) calculate local feature descriptor from interest points

% Get image gradient magnitudes and directions at each pixel

%  Convert gradient directions into corresponding bin from the 8 bins
%  (0:45:360)

% For each interest point, iterate over each cell in 3x3 grid (??? change
% to 32x32). For each cell, extract subarray of gradient magnitudes and
% directions. Calculate bin count weighted by gradient magnitudes for each
% cell. Concatente and normalize to get features for interest point

%  Return descriptors for each interest point



%% 
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