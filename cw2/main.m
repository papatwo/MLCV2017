%% Q1. Matching
close all; clear all
% 1) Manual
img = im2double(imread('C:\Users\akaishuu\Desktop\images.jpg'));
imshow(img);
if size(img)>2 % or selecting ONE colour channel
    img = rgb2gray(img);
end
blurMask=[0.03 0.105 0.222 0.286 0.222 0.105 0.03];
img = conv2(img, blurMask, 'same');
%[X1,Y1] = ginput(5); % get interest points in imgA
% [X2,Y2] = ginput(5); % get interest points in imgB

%% 2) Automatic
% a) Harris points detector
[dx, dy]= meshgrid(-1 :1,-1 : 1); % x & y derivative matrix
sigma=1;
% compute derivative of img (correlation between kernel and img)
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');
h = fspecial('gaussian',max(1, fix(6*sigma)),sigma); 
% autocorrelation matrix
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

%applying gaussian filter on the computed value (window thing)
Sx2 = conv2(Ix2, h, 'same');
Sy2 = conv2(Iy2, h, 'same');
Sxy = conv2(Ixy, h, 'same');
autoMat = [Sx2 Sxy; Sxy Sy2];

%%
% get the pixels col and row number
height = size(img,1);
width = size(img,2);
% set constant and threshold
k = 0.04; 
threshold = 0.007;

for i = 1:height
    for j = 1:width
        % Define at each pixel(i, j) the matrix M
        M = [Sx2(i,j) Sxy(i,j);Sxy(i,j) Sy2(i,j)];
        
        % Compute the response of the detector at each pixel (measure of
        % cornerness in terms of eve of M)
        % R<0: edge, R>0: corner, |R| small: flat
        R(i,j) = (det(M) - k * (trace(M)^2)); % equivalent to eve1*eve2-k(eve1+eve2)^2
    end
end
% or compute in once???
% R = (Sx2.*Sy2 - Sxy.^2)./(Sx2 + Sy2 + eps); % need a higher threshold

% find local max
radius = 1;
order = 2*radius +1;
mx = ordfilt2(R, 3^2,ones(3));
output = (R==mx)&(R>threshold);
[r, c]=find(output);
figure; imshow(img);hold on;plot(c,r,'ys') % c-xaxis; r-yaixs 
