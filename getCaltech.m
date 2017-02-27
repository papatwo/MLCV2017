function [ data_train, data_query ] = getCaltech( MODE )

showImg = 0; % Show training & testing images and their image feature vector (histogram representation)

PHOW_Sizes = [4 8 10]; % Multi-resolution, these values determine the scale of each layer.
PHOW_Step = 8; % The lower the denser. Select from {2,4,8,16}

switch MODE
    case 'Caltech' % Caltech dataset
        close all;
        imgSel = [15 15]; % randomly select 15 images each class without replacement. (For both training & testing)
        folderName = './Caltech_101/101_ObjectCategories';
        classList = dir(folderName);
        classList = {classList(3:end).name} % 10 classes
        
<<<<<<< HEAD
        numBins = 256; % for instance,
        
=======
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        disp('Loading training images...')
        % Load Images -> Description (Dense SIFT)
        cnt = 1;
        if showImg
            figure('Units','normalized','Position',[.05 .1 .4 .9]);
            suptitle('Training image samples');
        end
<<<<<<< HEAD
        
        desc_sel = cell(length(classList),1);
        codeword = cell(length(classList),1);
        
        figure;
        
=======
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        for c = 1:length(classList)
            subFolderName = fullfile(folderName,classList{c});
            imgList = dir(fullfile(subFolderName,'*.jpg'));
            imgIdx{c} = randperm(length(imgList));
            imgIdx_tr = imgIdx{c}(1:imgSel(1));
            imgIdx_te = imgIdx{c}(imgSel(1)+1:sum(imgSel));
            
            for i = 1:length(imgIdx_tr)
                I = imread(fullfile(subFolderName,imgList(imgIdx_tr(i)).name));
                
                % Visualise
<<<<<<< HEAD
                if i < 6 && showImg
=======
                if i < 6 & showImg
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
                    subaxis(length(classList),5,cnt,'SpacingVert',0,'MR',0);
                    imshow(I);
                    cnt = cnt+1;
                    drawnow;
                end
                
                if size(I,3) == 3
                    I = rgb2gray(I); % PHOW work on gray scale image
                end
                
                % For details of image description, see http://www.vlfeat.org/matlab/vl_phow.html
                [~, desc_tr{c,i}] = vl_phow(single(I),'Sizes',PHOW_Sizes,'Step',PHOW_Step); %  extracts PHOW features (multi-scaled Dense SIFT)
            end
<<<<<<< HEAD
            % Build visual vocabulary (codebook) for 'Bag-of-Words method'
            desc_sel{c} = single(vl_colsubset(cat(2,desc_tr{c,:}), 10e4)); % Randomly select 100k SIFT descriptors for clustering

            % K-means clustering
            disp(['Generating vocabulary for ' classList{c}]);
            [codeword{c}, visualword] = vl_kmeans(desc_sel{c}, numBins, 'Initialization', 'plusplus');
            % K++: greedily picks 10 maximally different data points as initial clustering centres
            
            subplot(5,2,c);
            histogram(visualword, numBins);
            title(classList{c});
            ylabel('Frequency');
            xlabel([classList{c} ' codeword']);        
        end
       
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,'\bf Visual Words Histograms for 10 image classes ','HorizontalAlignment', 'center','VerticalAlignment', 'top')
        
        
        % disp('Building visual codebook...')
=======
        end
        
        disp('Building visual codebook...')
        % Build visual vocabulary (codebook) for 'Bag-of-Words method'
        desc_sel = single(vl_colsubset(cat(2,desc_tr{:}), 10e4)); % Randomly select 100k SIFT descriptors for clustering
        
        % K-means clustering
        numBins = 256; % for instance,
        
        [centres, assignments] = vl_kmeans(desc_sel, 10, 'Initialization', 'plusplus');
        
        vocab = cell(10,1);
        for i = 1:10
            findWord = assignments == i;
            vocab{i} = desc_sel(:,findWord);
            [~,vocabSize] = size(vocab{i});
        end
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        
        
        
        
<<<<<<< HEAD
=======
       
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        disp('Encoding Images...')
        % Vector Quantisation
        
        % write your own codes here
        % ...
  
        
        % Clear unused varibles to save memory
        clearvars desc_tr desc_sel
end

switch MODE
    case 'Caltech'
        if showImg
<<<<<<< HEAD
            figure('Units','normalized','Position',[.05 .1 .4 .9]);
            suptitle('Testing image samples');
        end
        disp('Processing testing images...');
        cnt = 1;
        
        showSpatHist = 1; % Show visual word histogram of selected image
        
        codeword_te = cell(length(classList),length(imgIdx_te));
        
=======
        figure('Units','normalized','Position',[.05 .1 .4 .9]);
        suptitle('Testing image samples');
        end
        disp('Processing testing images...');
        cnt = 1;
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        % Load Images -> Description (Dense SIFT)
        for c = 1:length(classList)
            subFolderName = fullfile(folderName,classList{c});
            imgList = dir(fullfile(subFolderName,'*.jpg'));
            imgIdx_te = imgIdx{c}(imgSel(1)+1:sum(imgSel));
            
            for i = 1:length(imgIdx_te)
                I = imread(fullfile(subFolderName,imgList(imgIdx_te(i)).name));
                
                % Visualise
<<<<<<< HEAD
                if i < 6 && showImg
=======
                if i < 6 & showImg
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
                    subaxis(length(classList),5,cnt,'SpacingVert',0,'MR',0);
                    imshow(I);
                    cnt = cnt+1;
                    drawnow;
                end
                
                if size(I,3) == 3
                    I = rgb2gray(I);
                end
                [~, desc_te{c,i}] = vl_phow(single(I),'Sizes',PHOW_Sizes,'Step',PHOW_Step);
<<<<<<< HEAD
                [codeword_te{c,i}, visualword_te] = vl_kmeans(desc_te{c,i}, numBins, 'Initialization', 'plusplus');
                if showSpatHist
                    histogram(visualword_te, numBins);
                    title(['Visual Word Histogram for ' classList{c} ' Test Image']);
                    ylabel('Frequency');
                    xlabel('Test image codeword');                          
                end
            end
            

        end
        
        suptitle('Testing image samples');
        if showImg
            figure('Units','normalized','Position',[.5 .1 .4 .9]);
            suptitle('Testing image representations: 256-D histograms');
=======
            
            end
        end
        suptitle('Testing image samples');
                if showImg
            figure('Units','normalized','Position',[.5 .1 .4 .9]);
        suptitle('Testing image representations: 256-D histograms');
>>>>>>> 9f73a11a63a173c55b55e88c1a887905d1dbd26e
        end

        % Quantisation
        
        % write your own codes here
        % ...
        
        
    otherwise % Dense point for 2D toy data
        xrange = [-1.5 1.5];
        yrange = [-1.5 1.5];
        inc = 0.02;
        [x, y] = meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2));
        data_query = [x(:) y(:) zeros(length(x)^2,1)];
end
end

