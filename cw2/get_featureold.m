function features = get_featureold(img,intere_pt,patch_size)
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
r = intere_pt(2,:)';

length(c)
% As in coordinate convention: c-x, r-y
features = zeros(length(c),128); % 128 dimensions or 121???

% zero-pad image
img = im2double(imread(img));
% img = im2double(imread('test1.png'));
% img = im2double(imread('test2.png'));
if size(img)>2 % or selecting ONE colour channel
    img = rgb2gray(img);
end
blurMask=[0.03 0.105 0.222 0.286 0.222 0.105 0.03];
img = conv2(img, blurMask, 'same');
img_pad = padarray(img,[patch_size/2 patch_size/2],'symmetric');
[Gmag, Gdir] = imgradient(img_pad);

% Convert gradient direction into 8 bins
Gdir = Gdir + 180;
Gdir_bin = floor(Gdir/45)+1;

% iterate over each cell in 4x4 grid.
for i = 1:length(c) % the n-th interest point
    cnt = 0;
    for cellx_start = ceil(r(i) : patch_size / 4 : r(i)+patch_size -1) % each cell in 4x4 grid
        cellx_end = ceil(r(i)+patch_size/4-1);
        for celly_start = ceil(c(i) : patch_size / 4 : c(i) + patch_size - 1) % each cell in 4x4 grid
            celly_end = ceil(c(i)+patch_size/4-1);
            % extract magnitude and direction bin submatrices
            cell_Gmag = reshape(Gmag(celly_start:celly_end, cellx_start:cellx_end),[],1);
            cell_Gdir = reshape(Gdir_bin(celly_start:celly_end, cellx_start:cellx_end),[],1);
            % Compute features for each cell in grid
            for k = 1:length(cell_Gdir)
                idx = cnt*8+cell_Gdir(k); % feature idx (1-128)
                % Increments weigthed by gradient magnitude
                features(i, idx) = features(i, idx)+cell_Gmag(k); % prevent normalisation to be NaN
            end
            cnt=cnt+1;
        end         
    end
    features(i,:)=features(i,:)/norm(features(i,:));
end

end