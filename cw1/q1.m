init;
[data_train, data_test] = getData('Toy_Spiral'); 

%% to create  bags of dataset 
% 4 bags without replacement
% group = crossvalind('Kfold', length(data_train(:,1)), 4);
% for i=1:4
%     bag_uni{i} = data_train(group==i,:); 
% end

% bags with replacement
frac = 1 - 1/exp(1); % Bootstrap sampling fraction: 1 - 1/e (63.2%)
[N,D]=size(data_train);
bag={};
for i=1:4
    bag{i} = data_train(randsample(N,ceil(N*frac),1),:); % 4 bags of samples with replacement 
    figure; plot_toydata(bag{i}) 
end




%% grow trees 

%% Set the random forest parameters for instance, 
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try (what's this exactly mean???
param.split = 'IG';     % Currently support 'information gain' only

% grow all trees
trees = growTrees(data_train,param);
% prob field in trees struct: only the prob of first tree was calculated.
% to get all ten trees prob go to growTrees.m to line 47 and 49 to modify
% index of trees.....

%% Use Karpathy RF library
opts= struct;
opts.depth= 5;
opts.numTrees= 1;
opts.numSplits= 3;
opts.verbose= false;
opts.classifierID= 1; % weak learners to use. Can be an array for mix of weak learners too

X= bsxfun(@rdivide, bsxfun(@minus, data_train(:,1:2), mean(data_train(:,1:2))), var(data_train(:,1:2)));
Y=data_train(:,3);
frac = 1 - 1/exp(1); % Bootstrap sampling fraction: 1 - 1/e (63.2%)
[N,D]=size(data_train);
idx = randsample(N,ceil(N*frac),1);
q1_subX = X(idx,:);
q1_subY = Y(idx);
% classify
tic;
m= forestTrain(X, Y, opts);
timetrain= toc;
tic;

% for i=1:4
%     Y = m.treeModels{1,i}.weakModels{1,1}.Y;
%     YR = m.treeModels{1,i}.weakModels{1,1}.YR;
%     YL = m.treeModels{1,i}.weakModels{1,1}.YL;
%     
%     figure(i);
%     subplot(2,2,1)
%     bar(hist(Y,unique(Y)));
%     subplot(2,2,2)
%     bar(hist(YR,unique(YR)));
%     subplot(2,2,3)
%     bar(hist(YL,unique(YL)));
% end

yhatTrain = forestTest(m, X);
timetest= toc;

%%%change input to subbags data




%% QUESTION ONE FINAL
% axis-aligned weaker learner
load('q1_subbag.mat')
q1_subX= bsxfun(@rdivide, bsxfun(@minus, q1_subX(:,1:2), mean(q1_subX(:,1:2))), var(q1_subX(:,1:2)));
opts= struct;
opts.depth= 5;
opts.numTrees= 10;
opts.numSplits= 3;
opts.verbose= false;
ID=[1]; % choose different split function
%leafidx=randperm(15,4); %% WRONG!!! select 4 leafnodes to visualised their class distribution

leafidx=randperm(16,9); % can only select node 8-15 which are last internal nodes
% dist of some leafnodes are 0??? why????
data = [q1_subX q1_subY];
for i=1:length(ID) % loop different split funcs
tic;
opts.classifierID = ID(i);
tree{i}= forestTrain(data(:,1:2), data(:,end), opts); % get the tree for current split func
timetrain= toc;
tic;  
Y = tree{i}.treeModels{1,1}.weakModels{1,1}.Y; % parent node
YR = tree{i}.treeModels{1,1}.weakModels{1,1}.YR; % right children node
YL = tree{i}.treeModels{1,1}.weakModels{1,1}.YL; % left children node
dim = tree{i}.treeModels{1,1}.weakModels{1,1}.r; % selected axis
dec = tree{i}.treeModels{1,1}.weakModels{1,1}.dec; % threshold status
ig_best = tree{i}.treeModels{1,1}.weakModels{1,1}.maxgain;
t = tree{i}.treeModels{1,1}.weakModels{1,1}.t;
%XL = q1_subX(dec,:);
%XR = q1_subX(~dec,:);
leafdist = tree{i}.treeModels{1,1}.leafdist;

leaf=[];
for c=1:length(leafdist)
    [~,leaf_class] =max(leafdist(c,:)); % leaf nodes (actual counts of each model)
    leaf = [leaf;leaf_class]; 
end
% leafC{i} = leaf;
% figure;bar(hist(leafC{i},[1 2 3]))
% title(['Leaf distribution of SF' num2str(i) ' of all 16 leaf nodes of SF '])

figure; % class distribution of parent node and children nodes
% suptitle('Nodes Class Distribution')
subplot(2,2,1); 
bar(hist(Y,unique(Y)));title ('Class Distribution of Parent Node')
subplot(2,2,2)
bar(hist(YR,unique(Y)));title ('Right children Node')
subplot(2,2,3)
bar(hist(YL,unique(Y)));title ('Left children Node')
subplot(2,2,4)
RFvisualise_splitfunc(dec,data,dim,t,ig_best)
end
%%
% information gain of base node for tree constructed by different split
% funcs
for k=1:length(tree)
    Gain(k)=tree{1,k}.treeModels{1,1}.weakModels{1,1}.maxgain;
end

% class distribution of slected leaf nodes
figure;
for j=1:length(leafidx)
    subplot(3,3,j)
    bar(leafdist(j,:))
    title(sprintf('Leaf node %d',leafidx(j)))   
    %bar([1 2 3],leaf(leafidx(j),:));
    fprintf('class distribution of node %d. \n',leafidx(j))
end

%% GROW ALL TREES
%init;
%[data_train, data_test] = getData('Toy_Spiral'); 
X= bsxfun(@rdivide, bsxfun(@minus, data_train(:,1:2), mean(data_train(:,1:2))), var(data_train(:,1:2)));
Y= data_train(:,end);
opts= struct;
opts.depth= 5;
opts.numTrees= 10;
opts.numSplits= 3;
opts.verbose= false;
ID=[1]; % choose different split function
tree{1}= forestTrain(X, Y, opts); % get the tree for current split func


%% QUESTION TWO
% Test 4 data points given in the script
test_point = [-0.5 -0.7; 0.4 0.3; -0.7 0.4; 0.5 -0.5];

for n=1:length(test_point)%4 % how to get all 22801 test point prob
    leaves = testTrees([test_point(n,:) 0],trees);
    % average the class distributions of leaf nodes of all trees
    p_rf = trees(1).prob(leaves,:)
    p_rf_sum = sum(p_rf)/length(trees)
end
%%
for k=1:length(tree)
    [Yhard, Ysoft] = forestTest(tree{k}, test_point, opts);
    % Yhard are hard assignments to X's, Ysoft is NxK array of
    % probabilities, where there are K classes.
    % need to modify Ysoft if want prob as for first half question I changed
    % the prob to actual counts (in weakTrain)
    % Look at classifier distribution for fun, to see what classifiers were
    % chosen at split nodes and how often
    Yforwl{k}=[Yhard, Ysoft];
    fprintf('Classifier distributions %d:\n',k);
    classifierDist= zeros(1, 4);
    unused= 0;
    for i=1:length(tree{k}.treeModels)
        for j=1:length(tree{k}.treeModels{i}.weakModels)
            cc= tree{k}.treeModels{i}.weakModels{j}.classifierID;
            if cc>1 %otherwise no classifier was used at that node
                classifierDist(cc)= classifierDist(cc) + 1;
            else
                unused= unused+1;
            end
        end
    end
    fprintf('%d nodes were empty and had no classifier.\n', unused);
    for i=1:4
        fprintf('Classifier with id=%d was used at %d nodes.\n', i, classifierDist(i));
    end
    
    
    xrange = [-1.5 1.5];
    yrange = [-1.5 1.5];
    inc = 0.02;
    [x, y] = meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2));
    image_size = size(x);
    xy = [x(:) y(:)];
    
    [yhat, ysoft] = forestTest(tree{k}, xy);
    decmap= reshape(ysoft, [image_size 3]);
    decmaphard= reshape(yhat, image_size);
    
    subplot(121);
    imagesc(xrange,yrange,decmaphard);
    hold on;
    set(gca,'ydir','normal');
    cmap = [1 0.8 0.8; 0.95 1 0.95; 0.9 0.9 1];
    colormap(cmap);
    %%%%%%%%%%%%%%%modify here!=!!!!!!!!!!!!!!!
    plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
    plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
    plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
    hold off;
    title(sprintf('%d trees, Train time: %.2fs, Test time: %.2fs\n', opts.numTrees, timetrain, timetest));
    
    subplot(122);
    imagesc(xrange,yrange,decmap);
    hold on;
    set(gca,'ydir','normal');
    plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
    plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
    plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
    hold off;
    
    title(sprintf('Train accuracy: %f\n', mean(yhatTrain==Y)));
end






%%
subplot(2,2,4);
r = [-1.5 1.5]; % Data range
t = m.treeModels{1,1}.weakModels{1,1}.t;
if m.treeModels{1,1}.weakModels{1,1}.r == 1
    plot([t t],[r(1),r(2)],'r');
else
    plot([r(1),r(2)],[t t],'r');
end

hold on;
plot(data(data(:,end)==1,1), data(data(:,end)==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
hold on;
plot(data(data(:,end)==2,1), data(data(:,end)==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
hold on;
plot(data(data(:,end)==3,1), data(data(:,end)==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');

axis([r(1) r(2) r(1) r(2)]);
hold off;
%%

% Look at classifier distribution for fun, to see what classifiers were
% chosen at split nodes and how often
fprintf('Classifier distributions:\n');
classifierDist= zeros(1, 4);
unused= 0;
for i=1:length(m.treeModels)
    for j=1:length(m.treeModels{i}.weakModels)
        cc= m.treeModels{i}.weakModels{j}.classifierID;
        if cc>1 %otherwise no classifier was used at that node
            classifierDist(cc)= classifierDist(cc) + 1;
        else
            unused= unused+1;
        end
    end
end
fprintf('%d nodes were empty and had no classifier.\n', unused);
for i=1:4
    fprintf('Classifier with id=%d was used at %d nodes.\n', i, classifierDist(i));
end

%% plot results
xrange = [-1.5 1.5];
yrange = [-1.5 1.5];
inc = 0.02;
[x, y] = meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2));
image_size = size(x);
xy = [x(:) y(:)];

[yhat, ysoft] = forestTest(m, xy);
decmap= reshape(ysoft, [image_size 3]);
decmaphard= reshape(yhat, image_size);

subplot(121);
imagesc(xrange,yrange,decmaphard);
hold on;
set(gca,'ydir','normal');
cmap = [1 0.8 0.8; 0.95 1 0.95; 0.9 0.9 1];
colormap(cmap);
plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
hold off;
title(sprintf('%d trees, Train time: %.2fs, Test time: %.2fs\n', opts.numTrees, timetrain, timetest));

subplot(122);
imagesc(xrange,yrange,decmap);
hold on;
set(gca,'ydir','normal');
plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
hold off;