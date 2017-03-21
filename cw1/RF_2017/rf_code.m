function [ data_train,data_query, booktime, codebook ] = rf_code( MODE )
PHOW_Sizes = [4 8 10]; % Multi-resolution, these values determine the scale of each layer.
PHOW_Step = 8; % The lower the denser. Select from {2,4,8,16}

if strcmp(MODE, 'Caltech') % Caltech dataset
    close all;
    imgSel = [15 15]; % randomly select 15 images each class without replacement. (For both training & testing)
    folderName = './Caltech_101/101_ObjectCategories';
    classList = dir(folderName);
    classList = {classList(3:end).name} % 10 classes
    
    disp('Loading training images...')
    % Load Images -> Description (Dense SIFT)
    cnt = 1;
    
    for c = 1:length(classList)
        subFolderName = fullfile(folderName,classList{c});
        imgList = dir(fullfile(subFolderName,'*.jpg'));
        imgIdx{c} = randperm(length(imgList));
        imgIdx_tr = imgIdx{c}(1:imgSel(1));
        imgIdx_te = imgIdx{c}(imgSel(1)+1:sum(imgSel));
        
        for i = 1:length(imgIdx_tr)
            I = imread(fullfile(subFolderName,imgList(imgIdx_tr(i)).name));
            
            % Visualise
            
            
            if size(I,3) == 3
                I = rgb2gray(I); % PHOW work on gray scale image
            end
            
            % For details of image description, see http://www.vlfeat.org/matlab/vl_phow.html
            % each cell includes the descriptors of the image in that
            % class (collection of visual words is "Bag of 'words'")
            % desc: local area of an image (visual words)
            [~, desc_tr{c,i}] = vl_phow(single(I),'Sizes',PHOW_Sizes,'Step',PHOW_Step); %  extracts PHOW features (multi-scaled Dense SIFT)
            
        end
    end
    [ num_class, num_img ] = size(desc_tr);
    disp('Building visual codebook...')
    % Build visual vocabulary (codebook)
    % K-means clustering (codebook) for 'Bag-of-Words method'
    % codebook: independent features (a representative of similar patches)
    %         desc_sel = single(vl_colsubset(cat(2,desc_tr{:}), 10e4)); % Randomly select 100k SIFT descriptors for clustering
    desc_sel = [];
    for i = 1:size(desc_tr,1) % for all 10 classes
        desc_label = cat(2,desc_tr{i,:})'; % extract the training SIFT in current class
        idx = randsample(length(desc_label),10000); % randomly select 10k SIFT from this class
        desc_ = [desc_label(idx,:) i*ones(10000,1)]; % add class labels to the selected SIFT
        desc_sel = [desc_sel ; desc_]; % concatenate all selected classes together
    end
    tic
    param.num = 8;
    param.depth = 5;        % trees depth
    param.splitNum = 5;     % Number of split functions to try
    param.split = 'IG';
    desc_sel_tr = single(desc_sel);
    codebook = growTrees(desc_sel_tr,param);
    booktime=toc;
    
    % Vector Quantisation -- get data_train
    
    data_train = zeros(num_img*num_class,length(codebook(1).prob)+1);
    for i = 1:num_class
        for k=1:num_img
            dense_leaves = testTrees_fast(desc_tr{i,k}',codebook);
            dense_leaves(dense_leaves==0)=1;
            for c=1:length(codebook(1).prob) % num of nodes in this tree
                data_train(15*(i-1)+k,c)= length(find(dense_leaves==c));
            end
            data_train(15*(i-1)+k,end) = i;
        end
    end
    % Clear unused varibles to save memory
    clearvars desc_tr desc_sel
    
    
    %% processing testing images!!
    disp('Processing testing images...');
    cnt = 1;
    % Load Images -> Description (Dense SIFT)
    for c = 1:length(classList)
        subFolderName = fullfile(folderName,classList{c});
        imgList = dir(fullfile(subFolderName,'*.jpg'));
        imgIdx_te = imgIdx{c}(imgSel(1)+1:sum(imgSel));
        
        for i = 1:length(imgIdx_te)
            I = imread(fullfile(subFolderName,imgList(imgIdx_te(i)).name));
            
            % Visualise
            if size(I,3) == 3
                I = rgb2gray(I);
            end
            [~, desc_te{c,i}] = vl_phow(single(I),'Sizes',PHOW_Sizes,'Step',PHOW_Step);
            
        end
    end
    %         suptitle('Testing image samples');
    
    
    % Quantisation -- get data_test
    [ num_class, num_img ] = size(desc_te);
    data_query = zeros(num_class*num_img, length(codebook(1).prob)+1);
    for i = 1:num_class
        for k=1:num_img
            dense_leaves = testTrees_fast(desc_te{i,k}',codebook);
            for c=1:length(codebook(1).prob) % num of nodes in this tree
                data_query(15*(i-1)+k,c)= length(find(dense_leaves==c));
            end
            data_query(15*(i-1)+k,end) = i;
        end
    end
end
