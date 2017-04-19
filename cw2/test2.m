imgA = 'HG1.JPG';
imgB = 'HG2.JPG';
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
n = 4;
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

showMatchedFeatures(imgB, imgB, b, proj_ptB,'montage');
%%
abs(proj_ptB - testPoints2)
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
