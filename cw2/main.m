%% Q1. Matching

% 1) Manual
img = im2double(imread('C:\Users\akaishuu\Desktop\ÏÂÔØk.jpg'));
imshow(img);
if size(img)>2 % or selecting ONE colour channel
    img = rgb2gray(img);
end
%[X1,Y1] = ginput(5); % get interest points in imgA
% [X2,Y2] = ginput(5); % get interest points in imgB

%% 2) Automatic
% a) Harris points detector
dx = [-1 0 1;-1 0 1;-1 0 1]; % x derivative matrix
dy = dx'; % y derivative matrix

% compute derivative of img (correlation between kernel and img)
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');

% autocorrelation matrix
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

% apply Gaussian filter to autoMat
%applying gaussian filter on the computed value
h = fspecial('gaussian'); 
% autoMatConv = filter2(h,autoMat);
Ix2 = conv2(h,Ix2);
Iy2 = conv2(h,Iy2);
Ixy = conv2(h,Ixy);
autoMat = [Ix2 Ixy; Ixy Iy2];
% get the pixels col and row number
height = size(img,1);
width = size(img,2);
result = zeros(height,width); 
R = zeros(height,width);

alpha = 0.04;
Rmax = 0;
for i = 1:height
    for j = 1:width       
        M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)];
        R(i,j) = det(M) - alpha * (trace(M)^2);
        if R(i,j) > Rmax
            Rmax = R(i,j);
        end
    end
end

cnt = 0;
for i = 2:height-1
    for j = 2:width-1
        if R(i,j) > 0.1*Rmax && R(i,j) > R(i-1,j-1) && R(i,j) > R(i-1,j) && R(i,j) > R(i-1,j+1) && R(i,j) > R(i,j-1) && R(i,j) > R(i,j+1) && R(i,j) > R(i+1,j-1) && R(i,j) > R(i+1,j) && R(i,j) > R(i+1,j+1)
            result(i,j) = 1;
            cnt = cnt+1;
        end;
    end;
end;
[posc, posr] = find(result == 1);
cnt ;
imshow(img);
hold on;
plot(posr,posc,'r.');
