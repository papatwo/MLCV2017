%% Q1. Matching
% close all; clear all
% 1) Manual
img = im2double(imread('img1.pgm'));
imshow(img);
[X1,Y1] = ginput(8); % get interest points in imgA
[X2,Y2] = ginput(8); % get interest points in imgB
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
%%
imgA = 'img1.pgm';
imgB = 'img3.pgm';

imgA = im2double(imread(imgA));
if size(size(imgA),2)>2 % or selecting ONE colour channel
    imgA = rgb2gray(imgA);
end

imgB = im2double(imread(imgB));
if size(size(imgB),2)>2 % or selecting ONE colour channel
    imgB = rgb2gray(imgB);
end
%% Use self-written function
patch_size = 32;
Rthresh = 3000;
threshold = 0.95;

ptA = get_interePt(imgA, patch_size, Rthresh);
ptB = get_interePt(imgB, patch_size, Rthresh);

%featuresA = get_feature(imgA, ptA, patch_size);
%featuresB = get_feature(imgB, ptB, patch_size);
featuresA = get_features(imgA, ptA(1,:), ptA(2,:), patch_size);
featuresB = get_features(imgB, ptB(1,:), ptB(2,:), patch_size);

[matchmy, confidence,dist,r] = knn_match(featuresA, featuresB, threshold);

a = ptA(:,matchmy(:,1))';
b = ptB(:,matchmy(:,2))';

figure
showMatchedFeatures(imgA, imgB , a, b, 'montage');


%% Use builtin function
cornersA = detectHarrisFeatures(imgA);
figure;imshow(imgA); hold on;plot(cornersA.selectStrongest(50));

cornersB = detectHarrisFeatures(imgB);
figure;imshow(imgB); hold on;plot(cornersB.selectStrongest(20));

[features1,valid_points1] = extractFeatures(imgA,cornersA);
[features2,valid_points2] = extractFeatures(imgB,cornersB);
indexPairs = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
figure; showMatchedFeatures(imgA,imgB,matchedPoints1,matchedPoints2, 'montage');

%% 3) Transformation estimation

% SELF-WRITTEN: by using SVD solution to compute H (with mistakes)
H = get_homography(a,b);

% Implemented library code: (correct!!!)
a=ptA(:,matchmy(:,1))';
a=[a(1:4,:)';ones(1,4)];
b=ptB(:,matchmy(:,2))';
b=[b(1:4,:)';ones(1,4)];
H = homography2d(a, b);

% for transformation test:
i = 1; % change this number to any point you want to test
corre = H*[corre1(:,i);1];
corre = corre./corre(end); % norm the points vector
%% 
% b) fundamental matrix
load stereoPointPairs % for test convenience

%   input points for get H should be 2xN
[F,e1,e2] = get_fMatrix(matchedPoints1,matchedPoints2);
% % by using builtin
% F1 = estimateFundamentalMatrix(matchedPoints1,matchedPoints2);

%%
% c) projecting from imgB to A via homography.
%   calculate avg dist in pixels between these two mapped points
%   input points for get H should be Nx3
if size(matchedPoints1,2)<3
    H = homography2d([matchedPoints1(1:4,:) ones(4,1)]', [matchedPoints2(1:4,:) ones(4,1)]');
end

% for testing accuracy, find another set of match points
proj_ptB = [];
testPoints1=2%!!!;
testPoints2=2%!!!;
for i = 1:size(testPoints1,1)
    ptA = [testPoints1(i,:) 1];
    proj = H * ptA';
%     proj = proj./proj(end);
    proj(end)=[];
    proj = proj';
    proj_ptB = [proj_ptB; proj];
end
HA = sum(proj_ptB - testPoints2);

% MAYBE another way test
h11=H(1,1);h12=H(1,2);h13=H(1,3);
h31=H(3,1);h32=H(3,2);
x2pr=(h11*matchedPoints1(1,1)+h12*matchedPoints1(1,2)+h13)...
    /(h31*matchedPoints1(1,1)+h32*matchedPoints1(1,2)+1);
HA = sum(proj_ptB - matchedPoints2);

%%
% d) calculate epipolar line given point coor and a F

% Step through each matched pair of points and display the
% corresponding epipolar lines on the two images.

% x2'*F*x1 = 0 by fixing xy coor in x1 vector
% x2'* 3x1constant = 0;
const = Ft*[X1(1);Y1(1);1]; % epiline eq in img2 ?????

% calculate Epipolar line from F for every point in the image
% for test
load stereoPointPairs;
pt_1 = [matchedPoints1 ones(size(matchedPoints1,1),1)]';
pt_2 = [matchedPoints2 ones(size(matchedPoints2,1),1)]';

ll=[];lr=[];
figure;scatter(matchedPoints1(:,1),matchedPoints1(:,2));
figure;scatter(matchedPoints1(:,1),matchedPoints1(:,2));
for i = 1:size(pt_1,2)
    l2 = F*pt_1(:,i);    % Epipolar lines in image1: projection from imgA to B]'
    l1 = F'*pt_2(:,i);   % Epipolar lines in image2: projection from imgB to A
    points1=lineToBorderPoints(l2',[350 350]);  
    figure(2);hold on;line(points1(:,[1,3])',points1(:,[2,4])')
    points2=lineToBorderPoints(l1',[350 350]);
    figure(1);hold on;line(points2(:,[1,3])',points2(:,[2,4])')
end


% epline=epipolarLine(F',[X1 Y1])
% figure(1);hold on;line(points(:,[1,3])',points(:,[2,4])')
% epline2=epipolarLine(F',[X2 Y2])
% points2=lineToBorderPoints(epline,size(img2))
% figure(2);hold on;line(points2(:,[1,3])',points2(:,[2,4])')

%%
% Q2 Image Geometry
% 1) Homography
% a) HA before and after scale down the size of img
imgA = 'img1.pgm';
imgB = 'img2.pgm';

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

% Before Scale
test1 = cornersA.selectStrongest(50);
test2 = cornersB.selectStrongest(50);
test1 = test1.Location;
test2 = test2.Location;
if size(test1,2)<3
    H = homography2d([test1(1:4,:) ones(4,1)]', [test2(1:4,:) ones(4,1)]');
end

% for testing accuracy, find another set of match points
proj_ptB = [];
% testPoints1=!!!;
% testPoints2=!!!;
for i = 1:size(test1,1)
    ptA = [test1(i,:) 1];
    proj = H * ptA';
    proj = proj./proj(end);
    proj(end)=[];
    proj = proj';
    proj_ptB = [proj_ptB; proj];
end
HA = sum((proj_ptB - test2)./test2); % ?????how to normalise the error???????

% Scale down the size of img by 2
small_A = imresize (I1,0.5);
small_B = imresize (I2,0.5);
% detect scaled-down img intere_points
Spoints1 = detectHarrisFeatures(small_A);
Spoints2 = detectHarrisFeatures(small_B);
figure;imshow(small_A); hold on;plot(Spoints1.selectStrongest(50));
figure;imshow(small_B); hold on;plot(Spoints2.selectStrongest(50));
% pick up the first 50 strongest intere_pt
Stest1 = Spoints1.selectStrongest(50);
Stest2 = Spoints2.selectStrongest(50);
Stest1 = Stest1.Location;
Stest2 = Stest2.Location;
% compute new H matrix for scale-down img
if size(Stest1,2)<3
    H_s = homography2d([Stest1(1:4,:) ones(4,1)]', [Stest2(1:4,:) ones(4,1)]');
end

% for testing accuracy, find another set of match points
Sproj_ptB = [];
% testPoints1=!!!;
% testPoints2=!!!;
for i = 1:size(Stest1,1)
    ptA = [Stest1(i,:) 1];
    proj = H_s * ptA';
    proj = proj./proj(end);
    proj(end)=[];
    proj = proj';
    Sproj_ptB = [Sproj_ptB; proj];
end
HA_s = sum(Sproj_ptB - Stest2);

% img = im2double(imread('img1.pgm'));
% imshow(img);hold on;plot(X1,Y1,'ys');
% [X1,Y1] = ginput(8);
% img2 = im2double(imread('img2.pgm'));
% figure;imshow(img2);hold on;plot(X2,Y2,'ys');

%% 
% b) estimate H matrix from manually and automatically established corr
% better to clear var before run this section
% get correspondences manually
img = im2double(imread('img1.pgm'));
imshow(img);
[X1,Y1] = ginput(8);
hold on;plot(X1,Y1,'ys');

figure;
img2 = im2double(imread('img2.pgm'));
imshow(img2);
[X2,Y2] = ginput(8);
hold on;plot(X2,Y2,'ys');

x1 = [X1';Y1';ones(1,size(X1,1))];
% x1 = x1(:,1:4); % use only 4 points to compute H???
x2 = [X2';Y2';ones(1,size(X2,1))];
% x2 = x2(:,1:4);
manuH = homography2d(x1, x2)

% get correspondences automatically (currently using builtin)
imgA = 'img1.pgm';
imgB = 'img2.pgm';

I1 = im2double(imread(imgA));
I2 = im2double(imread(imgB));
points1 = detectHarrisFeatures(I1); % get intere_pt
points2 = detectHarrisFeatures(I2);
[features1,valid_points1] = extractFeatures(I1,points1); % get SIFT
[features2,valid_points2] = extractFeatures(I2,points2);
indexPairs = matchFeatures(features1,features2); % get correspondences
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

ax1 = [matchedPoints1.Location';ones(1,size(matchedPoints1.Location,1))];
ax1 = ax1(:,1:4); % use only 4 points to compute H???
ax2 = [matchedPoints2.Location';ones(1,size(matchedPoints2.Location,1))];
ax2 = ax2(:,1:4);
autoH = homography2d(ax1, ax2)

% visualise by existing code written by Peter Kovesi  
fig1=vgg_gui_H(imgA,imgB,manuH)
fig2=vgg_gui_H(imgA,imgB,autoH)

% analyse the geometirc transformation param!!!!!
%%%%%%%%%
%%%%%%%%%
%%%%%%%%%

%% 
% c) estimate H matrix from different no. correspondences (autoselect)
ax1all = [matchedPoints1.Location';ones(1,size(matchedPoints1.Location,1))];
ax2all = [matchedPoints2.Location';ones(1,size(matchedPoints2.Location,1))];
k=1;
for i = 4:30:length(ax1all)
    temp1 = ax1all(:,1:i); %
    temp2 = ax2all(:,1:i);
    Hdiff{k} = homography2d(temp1, temp2);
    k = k+1;
end

% compare HA and find outliers!!!!!!!!
%%%%%
%%%%%
%%%%%
%%%%%

%%
% 2) Stereo Vision (use images FD)
% a) estimate F matrix (using Q1.2a points)
% b) calculate epipoles and epipolar lines
ax1 = [matchedPoints1.Location';ones(1,size(matchedPoints1.Location,1))];
ax1 = ax1(:,1:8); % use only 4 points to compute H???
ax2 = [matchedPoints2.Location';ones(1,size(matchedPoints2.Location,1))];
ax2 = ax2(:,1:8);

[F,e1,e2] = get_fMatrix(ax1,ax2);
pt_1 = matchedPoints1.selectStrongest(8);
pt_2 = matchedPoints2.selectStrongest(8);
figure;imshow(imgA);hold on; plot(pt_1);
figure;imshow(imgB);hold on; plot(pt_2);
pt_1 = [(pt_1.Location)';ones(1,size(pt_1,1))];
pt_2 = [(pt_2.Location)';ones(1,size(pt_2,1))];

ll=[];lr=[];
for i = 1:size(pt_1,2)
    l2 = F*pt_1(:,i);    % Epipolar lines in image1: projection from imgA to B]'
    l1 = F'*pt_2(:,i);   % Epipolar lines in image2: projection from imgB to A
    points1=lineToBorderPoints(l2',[680 850]);  
    figure(2);hold on;line(points1(:,[1,3])',points1(:,[2,4])')
    points2=lineToBorderPoints(l1',[680 850]);
    figure(1);hold on;line(points2(:,[1,3])',points2(:,[2,4])')
end


%%
% c) calculate disparity map
imgA = imread('scene1.row3.col1.ppm');
imgB = imread('scene1.row3.col2.ppm');
disparityRange = [-6 10];
disparityMap = disparity(rgb2gray(imgA),rgb2gray(imgB),'BlockSize',...
    15,'DisparityRange',disparityRange);
figure
imshow(disparityMap,disparityRange);
title('Disparity Map');
colormap jet
colorbar
%%
% d) calculate and display depth map
disparityMap(abs(disparityMap)>100)=0;
figure;surf(disparityMap',imgA,'FaceColor','texturemap','EdgeColor','none');


%%
% e) change focal length by 2mm and do depth map again
%%%%%%
%%%%%%
%%%%%%

% add Gaussian noise to disparity map and do depth map again
Gaussian = randn(size(disparityMap));
Gaussian(abs(Gaussian)>2)=0;
disparity_noise = disparityMap+Gaussian;
figure;surf(disparity_noise',imgA,'FaceColor','texturemap','EdgeColor','none');
