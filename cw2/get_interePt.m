function intere_pt = get_interePt(img, patch_size,Rthresh)
% img: image path, must be string.
%% Q1. Matching
% close all; clear all
% 1) Manual
% img = im2double(imread('test1.png'));
% img = im2double(imread('test2.png'));
% imshow(img);

blurMask=[0.03 0.105 0.222 0.286 0.222 0.105 0.03];
% img = conv2(img, blurMask, 'same');
%[X1,Y1] = ginput(5); % get interest points in imgA
% [X2,Y2] = ginput(5); % get interest points in imgB

%% 2) Automatic
% a) Harris points detector
[dx, dy]= meshgrid(-1 :1,-1 : 1); % x & y derivative matrix

% compute derivative of img (correlation between kernel and img)
Ix = imfilter(img, dx);
Iy = imfilter(img, dy);
 
% autocorrelation matrix
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

%applying gaussian filter on the computed value (window thing)
% h = fspecial('gaussian',max(1, fix(6*sigma)),sigma);
h = blurMask;
Sx2 = imfilter(Ix2, h);
Sy2 = imfilter(Iy2, h);
Sxy = imfilter(Ixy, h);

% get the pixels col and row number
% height = size(img,1);
% width = size(img,2);
% set constant and threshold
% k = 0.04; 
% threshold = 0.007;
% for i = 1:height
%     for j = 1:width
%         % Define at each pixel(i, j) the matrix M
%         M = [Sx2(i,j) Sxy(i,j);Sxy(i,j) Sy2(i,j)];
%         
%         % Compute the response of the detector at each pixel (measure of
%         % cornerness in terms of eve of M)
%         % R<0: edge, R>0: corner, |R| small: flat
%         R(i,j) = (det(M) - k * (trace(M)^2)); % equivalent to eve1*eve2-k(eve1+eve2)^2
%     end
% end

alpha = 0.04;
R = (Sx2.*Sy2 - Sxy.^2)-alpha*(Sx2 + Sy2).^2;
R_sort = sort(reshape(R,[],1),'descend');
% <<<<<<< HEAD
% threshold = R_sort(50);
% threshold = 0.7;
% =======
threshold = mean(R_sort(1:min([Rthresh,length(R_sort)])));
% >>>>>>> acf1ea56d5e6892d834917b00b1393813224e5da
% Remove low gardients, graythresh makes the threshold adapt to image
highR = R>threshold;

% find local max: Non-maximum suppression, finds next highest point in
% neighborhood and checks whether this point in R is greater than that
nearest_max = ordfilt2(R, 8, [1 1 1;1 0 1;1 1 1]);
local_maxima = R>nearest_max;

% get interest points: satisfy both high gradient and is local maximum
% around the neighborhood
output = highR & local_maxima;

% radius = 1;
% order = 2*radius +1;
% mx = ordfilt2(R, 3^2,ones(3));
% output = (R>mx)&(R>threshold);

% remove interest points near image boundary
boundary = patch_size+1;
output(1:boundary,:)=0;
output(end-boundary:end,:)=0;
output(:,1:boundary)=0;
output(:,end-boundary:end)=0;

[r, c]=find(output);
intere_pt=[c';r']; % return interest points matrix
figure; imshow(img);
hold on;plot(c,r,'ys') % c-xaxis; r-yaixs find 69 interst points
end