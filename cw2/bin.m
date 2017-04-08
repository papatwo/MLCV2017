%% Q1. Matching
% close all; clear all
% 1) Manual
img = im2double(imread('images.jpg'));
% img = im2double(imread('test1.png'));
% img = im2double(imread('test2.png'));
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

% get the pixels col and row number
height = size(img,1);
width = size(img,2);
% set constant and threshold
k = 0.04; 
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

% % or compute in once???
R = (Sx2.*Sy2 - Sxy.^2)./(Sx2 + Sy2 + eps); % need a higher threshold
threshold = 0.1*graythresh(R);

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
boundary = 6*2;
output(1:boundary,:)=0;
output(end-boundary:end,:)=0;
output(:,1:boundary)=0;
output(:,end-boundary:end)=0;

[r, c]=find(output);
intere_pt=[c';r']; % interest points matrix
figure; imshow(img);
hold on;plot(c,r,'ys') % c-xaxis; r-yaixs find 69 interst points

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

% As in coordinate convention: c-x, r-y
features = zeros(length(c),128); % 128 dimensions or 121???

% zero-pad image
img_pad = padarray(img,[6/2 6/2],'symmetric');
[Gmag, Gdir] = imgradient(img_pad);

% Convert gradient direction into 8 bins
Gdir = Gdir + 180;
Gdir_bin = floor(Gdir/45)+1;

% iterate over each cell in 3x3 grid.
for i = 1:length(c) % the n-th interest point
    cnt = 0;
    for cellx_start = r(i) : 6 / 3 : r(i)+6 -1 % each cell in 4x4 grid
        cellx_end = r(i)+6/3-1;
        for celly_start = c(i) : 6 / 3 : c(i) + 6 - 1 % each cell in 4x4 grid
            celly_end = c(i)+6/3-1;
            % extract magnitude and direction bin submatrices
            cell_Gmag = reshape(Gmag(celly_start:celly_end, cellx_start:cellx_end),[],1);
            cell_Gdir = reshape(Gdir_bin(celly_start:celly_end, cellx_start:cellx_end),[],1);
            % Compute features for each cell in grid
            for k = 1:length(cell_Gdir)
                idx = cnt*8+cell_Gdir(k);
                % Increments weigthed by gradient magnitude
                features(i, idx) = features(i, idx)+cell_Gmag(k); % prevent normalisation to be NaN
            end
            cnt=cnt+1;
        end         
    end
    features(i,:)=features(i,:)/norm(features(i,:));
end
%%
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

