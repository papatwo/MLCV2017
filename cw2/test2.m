imgA = 'HG1.JPG';
imgB = 'HG2.JPG';
C2 = imgB;
% imgA = 'img1.pgm';
% imgB = 'img2.pgm';
imgA = im2double(imread(imgA));
if size(size(imgA),2)>2 % or selecting ONE colour channel
    imgA = rgb2gray(imgA);
end

imgB = im2double(imread(imgB));
if size(size(imgB),2)>2 % or selecting ONE colour channel
    imgB = rgb2gray(imgB);
end


patch_size = 32;
threshold = 0.87;
Rthresh = 3000;

% a = csvread('a.dat');
% b = csvread('b.dat');
% ptA = csvread('ptA.dat');
% ptB = csvread('ptB.dat');
% featuresA= csvread('featuresA.dat');
% featuresB = csvread('featuresB.dat');
% matchmy = csvread('matchmy.dat');
% 
% 
ptA = get_interePt(imgA, patch_size, Rthresh);
ptB = get_interePt(imgB, patch_size, Rthresh);

% featuresA = get_feature(I1,ptA, patch_size);
% featuresB = get_feature(I2, ptB, patch_size);
featuresA = get_features(imgA, ptA(1,:), ptA(2,:), patch_size);
featuresB = get_features(imgB, ptB(1,:), ptB(2,:), patch_size);
[matchmy, confidence,dist,r] = knn_match(featuresA, featuresB, threshold);

[~,L_sort] = sort(confidence, 2);
matchmy = matchmy(L_sort',:);

a = ptA(:,matchmy(:,1))';
b = ptB(:,matchmy(:,2))';
figure
showMatchedFeatures(imgA, imgB, a, b, 'montage');


%% 3) Transformation estimation

% % SELF-WRITTEN: by using SVD solution to compute H (with mistakes)
% %   input points for get H should be 2xN 
% % 4 matched pts are enough for the computation
% H = get_homography(a(1:4,:)',b(1:4,:)');
% 
% % % Implemented library code: (correct!!!)
% % a=[a(1:4,:)';ones(1,4)];
% % b=[b(1:4,:)';ones(1,4)];
% % H = homography2d(a, b);
% % H = H./H(end,end);
% 
% % for transformation test:
% i = 3; % change this number to any point you want to test
% corre = H*[a(i,:)';1];
% corre = corre./corre(end);

%%
% c) projecting from imgB to A via homography.
%   calculate avg dist in pixels between these two mapped points
%   input points for get H should be Nx3
% if size(matchedPoints1,2)<3
%     H = homography2d([matchedPoints1(1:4,:) ones(4,1)]', [matchedPoints2(1:4,:) ones(4,1)]');
% end
n = 8;
format long g;
% 
% imgA = 'HG1.JPG';
% imgB = 'HG3.JPG';
% imgA =rgb2gray(im2double(imread(imgA)));
% imgB =rgb2gray(im2double(imread(imgB)));
% 
% X1 = csvread('X1.dat');
% X2 = csvread('X3.dat');
% Y1 = csvread('Y1.dat');
% Y2 = csvread('Y3.dat');
% 
% a = [X1 Y1];
% b = [X2 Y2];

H = get_homography(a(1:n,:)',b(1:n,:)') 

proj_ptB = [];
%testPoints1=a(1:4,:);%!!!;
testPoints1=a;
testPoints2=b;
%testPoints2=b(1:4,:);%!!!;
for i = 1:size(testPoints1,1)
    ptA = [testPoints1(i,:) 1];
    proj = H * ptA';
    proj = proj./proj(end);
    proj(end)=[];
    proj = proj';
    proj_ptB = [proj_ptB; proj];
end
HA = sum(abs(proj_ptB - testPoints2))./size(a,1)
figure;
showMatchedFeatures(imgB, imgB, b, proj_ptB,'montage');


%% outlier!!!!!!!!!!!!!!!
% 
% idx=abs(proj_ptB - testPoints2)<10;
% ori_inlier=testPoints2(sum(idx,2)>0,:)
% proj_inlier=proj_ptB(sum(idx,2)>0,:)
% figure;
% showMatchedFeatures(imgB, imgB, proj_inlier, ori_inlier,'montage');
%% Outlier detection
% Compare histogram of hue of interest point areas in feature space
% Threshold distance to determine inlier/outlier

colourB = im2double(imread(C2));
[numR, numC] = size(colourB(:,:,1));
shift = patch_size/2; % Use patch size to determine area size
bin = 360; % Hue as 360 degrees
colourDist = zeros(size(testPoints2,1), 1);

for i = 1 : size(testPoints2,1)
    left_boundB = testPoints2(i,1)-shift;
    right_boundB = testPoints2(i,1)+shift;
    up_boundB = testPoints2(i,2)-shift;
    low_boundB = testPoints2(i,2)+shift;
    
    left_boundProj = proj_ptB(i,1)-shift;
    right_boundProj = proj_ptB(i,1)+shift;
    up_boundProj = proj_ptB(i,2)-shift;
    low_boundProj = proj_ptB(i,2)+shift;
    
    % reduce patch size if area exceed image
    if left_boundB < 1 || up_boundB < 1 || left_boundProj < 1 || up_boundProj < 1
        reduce = min([left_boundB, up_boundB, left_boundProj, up_boundProj]);
        left_boundB = left_boundB + reduce;
        right_boundB = right_boundB - reduce;
        up_boundB = up_boundB + reduce;
        low_boundB = low_boundB - reduce;

        left_boundProj = left_boundProj + reduce;
        right_boundProj = right_boundProj - reduce;
        up_boundProj = up_boundProj + reduce;
        low_boundProj = low_boundProj - reduce;
    end
    
    if right_boundB > numC || low_boundB > numR || right_boundProj  > numC || low_boundProj   > numR 
        reduce = max([right_boundB - numC, low_boundB - numR, right_boundProj - numC, low_boundProj- numR]);
        left_boundB = left_boundB + reduce;
        right_boundB = right_boundB - reduce;
        up_boundB = up_boundB + reduce;
        low_boundB = low_boundB - reduce;

        left_boundProj = left_boundProj + reduce;
        right_boundProj = right_boundProj - reduce;
        up_boundProj = up_boundProj + reduce;
        low_boundProj = low_boundProj - reduce;
    end
    
    
    colourPatchB = colourB(up_boundB : low_boundB, left_boundB : right_boundB, :);
    colourPatchProj = colourB(up_boundProj : low_boundProj, left_boundProj: right_boundProj, :);
    HSV_B = rgb2hsv(colourPatchB);
    HSV_proj = rgb2hsv(colourPatchProj);
    HueB = reshape(HSV_B(:,:,1),1,[]);
    HueProj = reshape(HSV_proj(:,:,1),1,[]);
    H_histB = histcounts(HueB * 360 , bin); % Histgram of hue
    H_histProj = histcounts(HueProj * 360 , bin); % Histgram of hue
    colourDist(i,1) = pdist2(H_histB, H_histProj);
end

C_thresh = 5;

outliers = find(colourDist > C_thresh);
inliers = find(colourDist < C_thresh);

figure;showMatchedFeatures(colourB, colourB, b(inliers,:), proj_ptB(inliers,:),'montage');
figure;showMatchedFeatures(colourB, colourB, b(outliers,:), proj_ptB(outliers,:),'montage');


%%
% % % MAYBE another way test
% h11=H(1,1);h12=H(1,2);h13=H(1,3);
% h31=H(3,1);h32=H(3,2);
% x2pr=(h11*matchedPoints1(1,1)+h12*matchedPoints1(1,2)+h13)...
%     /(h31*matchedPoints1(1,1)+h32*matchedPoints1(1,2)+1);
% HA = sum(proj_ptB - matchedPoints2);

%%
% Q2 Image Geometry
% 1) Homography
% a) HA before and after scale down the size of img
clear all;
clc;

imgA = 'HG1.JPG';

imgA = im2double(imread(imgA));
if size(size(imgA),2)>2 % or selecting ONE colour channel
    imgA = rgb2gray(imgA);
end

patch_size = 32;
Rthresh = 3000;
threshold = 0.87;
small_B = imresize (imgA,0.5);
ptA = get_interePt(imgA, patch_size, Rthresh);
ptB = get_interePt(small_B, patch_size, Rthresh);
featuresA = get_features(imgA, ptA(1,:), ptA(2,:), patch_size);
featuresB = get_features(small_B, ptB(1,:), ptB(2,:), patch_size);

[matchmy, confidence,dist,r] = knn_match(featuresA, featuresB, threshold);

[~,L_sort] = sort(confidence, 2);

matchmy = matchmy(L_sort',:);
a = ptA(:,matchmy(:,1))';
b = ptB(:,matchmy(:,2))';


n = 7;
showMatchedFeatures(imgA, small_B, a(1:n,:), b(1:n,:));

H = get_homography(a(1:n,:)',b(1:n,:)') 
 % when H is computed from 15 pts it has lowest projection error!!!!!!!!!!!
% for testing accuracy, find another set of match points1
% 
% figure; imshow(imgA);
% hold on;plot(a(1:n,1),a(1:n,2),'ys','Color', 'r');
% 
% figure; imshow(imgB);
% hold on;plot(b(1:n,1),b(1:n,2),'ys');

proj_ptB = [];
%testPoints1=a(1:4,:);%!!!;
testPoints1=aa;
testPoints2=ab;
%testPoints2=b(1:4,:);%!!!;
for i = 1:size(testPoints1,1)
    ptA = [testPoints1(i,:) 1];
    proj = autoH * ptA';
    proj = proj./proj(end);
    proj(end)=[];
    proj = proj';
    proj_ptB = [proj_ptB; proj];
end
HA = sum(abs(proj_ptB - testPoints2))./size(aa,1)



% img = im2double(imread('img1.pgm'));
% imshow(img);hold on;plot(X1,Y1,'ys');
% [X1,Y1] = ginput(8);
% img2 = im2double(imread('img2.pgm'));
% figure;imshow(img2);hold on;plot(X2,Y2,'ys');

%% 
% b) estimate H matrix from manually and automatically established corr
% better to clear var before run this section
% get correspondences manually
imgA = 'HG3.JPG';
imgB = 'HG1.JPG';

img = im2double(imread(imgA));
imshow(img);
[X1,Y1] = ginput(16);
hold on;plot(X1,Y1,'ys');

figure;
img2 = im2double(imread(imgB));
imshow(img2);
[X2,Y2] = ginput(16);
hold on;plot(X2,Y2,'ys');

n=8;
x1 = [X1';Y1'];
manux1 = x1(:,1:n); % use only 4 points to compute H???
x2 = [X2';Y2'];
manux2 = x2(:,1:n);
manuH = get_homography(manux1,manux2)

% get correspondences automatically (currently using builtin)
imgA = 'HG1.JPG';
imgB = 'HG2.JPG';

imgA = im2double(imread(imgA));
if size(size(imgA),2)>2 % or selecting ONE colour channel
    imgA = rgb2gray(imgA);
end

imgB = im2double(imread(imgB));
if size(size(imgB),2)>2 % or selecting ONE colour channel
    imgB = rgb2gray(imgB);
end
patch_size = 32;
Rthresh = 3000;
threshold = 0.87;
n = 8;

aptA = get_interePt(imgA, patch_size, Rthresh);
aptB = get_interePt(imgB, patch_size, Rthresh);
afeaturesA = get_features(imgA, aptA(1,:), aptA(2,:), patch_size);
afeaturesB = get_features(imgB, aptB(1,:), aptB(2,:), patch_size);

showMatchedFeatures(imgA, imgB, a(1:n,:), b(1:n,:));

[amatchmy, confidence,dist,r] = knn_match(afeaturesA, afeaturesB, threshold);

[~,L_sort] = sort(confidence, 2);
amatchmy = amatchmy(L_sort',:);

aa = aptA(:,amatchmy(:,1))';
ab = aptB(:,amatchmy(:,2))';
autoH = get_homography(aa(1:n,:)',ab(1:n,:)');

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
