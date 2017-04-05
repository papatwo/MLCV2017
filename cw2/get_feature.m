function features = get_feature(img,intere_pt,patch_size)
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
c = intere_pt(1,:)';
% r = intere_pt(2,:)';

length(c)
% As in coordinate convention: c-x, r-y
features = zeros(length(c),128); % 128 dimensions or 121???

% zero-pad image
img = im2double(imread(img));
% img = im2double(imread('test1.png'));
% img = im2double(imread('test2.png'));
% if size(size(img),2)>2 % or selecting ONE colour channel
%     img = rgb2gray(img);
% end
% blurMask=[0.03 0.105 0.222 0.286 0.222 0.105 0.03];
% img = conv2(img, blurMask, 'same');
img_pad = padarray(img,[patch_size/2 patch_size/2],'symmetric');
[Gmag, Gdir] = imgradient(img_pad,'central');

% Convert gradient direction into 8 bins
Gdir = Gdir + 180;
Gdir_bin = floor(Gdir/45)+1;

% iterate over each cell in 4x4 grid.
hf = patch_size/2;
x = intere_pt(1,:);
y = intere_pt(2,:);
for i = 1:length(c) % the n-th interest point
%     cnt = 0;
    
    Xi = x(i)-hf; % patch Xcoor start
    Xe = x(i)+hf-1; % patch Xcoor end
    Yi = y(i)-hf; % patch Ycoor start
    Ye = y(i)+hf-1; % patch Ycoor end
    %find all mag&dir of the 32x32 patch around intere_pt
    patch_mag = Gmag(Yi:Ye,Xi:Xe); % because row-ycoor & col
    patch_dirBin = Gdir_bin(Yi:Ye,Xi:Xe);
    
    cell_Gmag = im2col(patch_mag,[hf/2 hf/2],'distinct');
    cell_dirBin = im2col(patch_dirBin,[hf/2 hf/2],'distinct');
    vector = zeros(8, size(cell_dirBin,2)); % in each patch has 16 cells
    for t = 1:size(cell_Gmag,2) % loop through 16 cells in one patch
       [~,~,ind] = histcounts(cell_dirBin(:,t),8); 
       for s = 1:size(cell_Gmag,1) % loop through 64 points in one cell
        vector(ind(s),t) = features(ind(s),t)+cell_Gmag(s,t); % every count on a particular bin add the mag
       end
    end
    revector = reshape(vector,[1 128]); % expand the 128 features at this point
    revector = revector / norm(revector); 
    features(i,:) = revector; % concatenate with the summary feature matrix
end    
    
%     for k = 1:4
%         x_begin_grid = Xi+(k-1)*patch_size/4;
%         x_end_grid = Xi+k*patch_size/4;
%         for n = 1:4
%             y_begin_grid = Yi+(n-1)*patch_size/4;
%             y_end_gird = Yi+n*patch_size/4;
%             cell_Gmag = reshape(patch_mag(x_begin_grid:x_end_grid,y_begin_grid:y_end_gird),[],1);
%             cell_Gdir = reshape(patch_dirBin(x_begin_grid:x_end_grid,y_begin_grid:y_end_gird),[],1);
%             for c = 1:length(cell_Gdir)
%                 idx = cnt*8+cell_Gdir(c); % feature idx (1-128)
%                 % Increments weigthed by gradient magnitude
%                 features(i, idx) = features(i, idx)+cell_Gmag(c); % prevent normalisation to be NaN
%             end 
%             cnt = cnt+1;
%         end
%     end
% 
% %     for cellx_start = r(i)-patch_size/2 : patch_size / 4 : r(i)+patch_size/2  % each cell in 4x4 grid
% %         cellx_end = r(i)+patch_size/4-1;
% %         for celly_start = c(i)-patch_size/2 : patch_size / 4 : c(i) + patch_size/2 - 1 % each cell in 4x4 grid
% %             celly_end = c(i)+patch_size/4-1;
% %             % extract magnitude and direction bin submatrices
% %             cell_Gmag = reshape(Gmag(celly_start:celly_end, cellx_start:cellx_end),[],1);
% %             cell_Gdir = reshape(Gdir_bin(celly_start:celly_end, cellx_start:cellx_end),[],1);
% %             % Compute features for each cell in grid
% %             for k = 1:length(cell_Gdir)
% %                 idx = cnt*8+cell_Gdir(k); % feature idx (1-128)
% %                 % Increments weigthed by gradient magnitude
% %                 features(i, idx) = features(i, idx)+cell_Gmag(k); % prevent normalisation to be NaN
% %             end
% %             cnt=cnt+1;
% %         end         
% %     end
%     features(i,:)=features(i,:)/norm(features(i,:));
% end


