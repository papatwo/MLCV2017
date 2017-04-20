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
%%
testPoints1=a;
testPoints2=b;
% HA=zeros(length(4:6:35),2); 
k = 1;
HA=[]; proj_cell = {};
for n = 4:6:35    
    H = get_homography(a(1:n,:)',b(1:n,:)');
    proj_ptB = [];
    for i = 1:size(testPoints1,1)
        ptA = [testPoints1(i,:) 1];
        proj = H * ptA';
        proj = proj./proj(end);
        proj(end)=[];
        proj = proj';
        proj_ptB = [proj_ptB; proj];
    end
    HA(k,:) = sum(abs(proj_ptB - testPoints2))./size(a,1);
    proj_cell{k} = proj_ptB;
    k=k+1;
end

for aa = 1:6  
    subplot(2,3,aa); plot(testPoints2(:,1),testPoints2(:,2),'ro')
    hold on ; plot(proj_cell{aa}(:,1),proj_cell{aa}(:,2),'gs')   
end

legend()