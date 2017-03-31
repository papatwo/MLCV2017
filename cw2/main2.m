%% Q1. Matching

% 1) Manual
img = im2double(imread('C:\Users\akaishuu\Desktop\images.jpg'));
imshow(img);
img=img(:,:,1);
sigma=1;
%[X1,Y1] = ginput(5); % get interest points in imgA
% [X2,Y2] = ginput(5); % get interest points in imgB

%% 2) Automatic
% a) Harris points detector
[dx, dy] = meshgrid(-1 :1,-1:1); % x derivative matrix
% compute derivative of img (correlation between kernel and img)
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');

h = fspecial('gaussian',max(1, fix(6*sigma)),sigma); 
% autocorrelation matrix
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

Sx2 = conv2(Ix2, h, 'same');
Sy2 = conv2(Iy2, h, 'same');
Sxy = conv2(Ixy, h, 'same');
%%
height = size(img,1);
width = size(img,2);
for i = 1:height
    for j = 1:width 
        % Define at each pixel(i, j) the matrix M
        M = [Sx2(i,j) Sxy(i,j);Sxy(i,j) SIy2(i,j)];
        
        % Compute the response of the detector at each pixel (measure of
        % cornerness in terms of eve of M)
        % R<0: edge, R>0: corner, |R| small: flat
        R(i,j) = (det(M) - k * (trace(M).^2));
        
    end
end
R = (Sx2.*Sy2 - Sxy.^2)./(Sx2 + Sy2 + eps);
radius = 1;
order = 2*radius +1;
mx = ordfilt2(R, order^2,ones(order));
harris_points=(R==mx)&(R>0.07);
[r c]=find(harris_points)
figure; imshow(img);hold on;plot(c,r,'ys')


