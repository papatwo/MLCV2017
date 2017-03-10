% Simple Random Forest Toolbox for Matlab
% written by Mang Shao and Tae-Kyun Kim, June 20, 2014.
% updated by Tae-Kyun Kim, Feb 09, 2017

% This is a guideline script of simple-RF toolbox.
% The codes are made for educational purposes only.
% Some parts are inspired by Karpathy's RF Toolbox

% Under BSD Licence

%% Initialisation
init;

% Select dataset
[data_train, data_test] = getData('Toy_Spiral'); % {'Toy_Gaussian', 'Toy_Spiral', 'Toy_Circle', 'Caltech'}

% Set the random forest parameters for instance, 
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only



%%%%%%%%%%%%%
% check the training and testing data
    % data_train(:,1:2) : [num_data x dim] Training 2D vectors
    % data_train(:,3) : [num_data x 1] Labels of training data, {1,2,3}
    
plot_toydata(data_train);

    % data_test(:,1:2) : [num_data x dim] Testing 2D vectors, 2D points in the
    % uniform dense grid within the range of [-1.5, 1.5]
    % data_train(:,3) : N/A
    
scatter(data_test(:,1),data_test(:,2),'.b');


%% % Train Random Forest
%% Q1 create subbags
createSubbag;
%% Q2 split the FIRST node
axisAlignSplitFirstNode;

%% CHANGE SPLIT FUNCTION???

for i=1:10
% splitting the first node
bestig= -100;
iter = param.splitNum;
data = bag{1,1}; % use one subbag to split one tree
[N,D]=size(data);
for n = 1:iter
    visualise = 0;
    % axis-aligned split function:
    r1 = randi(D-1);
    r2 = randi(D-1); % Pick one random dimension -- pick x-axis or y-axis
    w= randn(3, 1); %
    idx_ = [data(:, [r1 r2]), ones(N, 1)]*w < 0; 
    dim = [r1 r2];
    t = w;
    l_data = data_train(idx_,:); % data go LEFT branch
    r_data = data_train(~idx_,:); % data go RIGHT branch
    
    ig = getIG(data_train, idx_);
    
    if bestig < ig
        bestig = ig;
        t_best = t;
        dim_best = dim;
        idx_best = idx_;
    end
    
%     if visualise
%         visualise_splitfunc(idx_,data_train,dim,t,bestig,0);
%         pause();
%     end    
end

figure; % visualise the best split
visualise_splitfuncTest(idx_best,data_train,dim_best,t_best,bestig,0);
end
%% Q3 Grow ONE complete tree
% for growing 1 tree first, reset the param
% T = 1;
% idx=1:length(bag{1}); 
% trees(T).node(1) = struct('idx',idx,'t',nan,'dim',-1,'prob',[]);
% 
% % Split Nodes
% for n = 1:2^(param.depth-1)-1
%     [trees(T).node(n),trees(T).node(n*2),trees(T).node(n*2+1)] = splitNode(data_train,trees(T).node(n),param);
% end

param.num=1;
trees=growTrees(data_train,param);
% Store class distribution in the leaf nodes.
makeLeaf;
% Visualise the leaf node class distribution
figure;
visualise_leaf
% every time train a tree will have different number of leaf node

%% linear grow one tree
% T = 1;
% idx=1:length(bag{1}); 
% trees(T).node(1) = struct('idx',idx,'t',nan,'dim',-1,'prob',[]);
% 
% % Split Nodes
% for n = 1:2^(param.depth-1)-1
%     [trees(T).node(n),trees(T).node(n*2),trees(T).node(n*2+1)] = splitNodeLinear(data_train,trees(T).node(n),param);
% end
param.num=1;
trees=growTreesLinear(data_train,param);
% Store class distribution in the leaf nodes.
makeLeaf;
% Visualise the leaf node class distribution
figure;
visualise_leaf

%% Q4 Grow all 10 trees
param.num=10;

trees_axis = growTrees(data_train,param); 

trees_linear=growTreesLinear(data_train,param);
% growTree is currently using linear split 3.9 16:41
%% Q5 Testing data points
%%%%%%%%%%%%%%%%%%%%%%
% Evaluate/Test Random Forest
pointTestTreeAxis;

%%
% grab the few data points and evaluate them one by one by the leant RF
test_point = [-.5 -.7; .4 .3; -.7 .4; .5 -.5];
figure;plot_toydata(data_train);
plot(test_point(:,1), test_point(:,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.9 .9 .9], 'MarkerEdgeColor','k');
title('Visualisation of Test Points');

for n=1:length(test_point)%4 % how to get all 22801 test point prob
    leaves = testTreesLinear([test_point(n,:) 0],trees_linear);
    % average the class distributions of leaf nodes of all trees
    p_rf = trees_linear(1).prob(leaves,:);
    p_rf_sum(n,:) = mean(p_rf,1);
    figure;
    for p = 1:length(leaves) % visualise the leaf class distribution of all ten trees
        subplot(2,5,p);bar(p_rf(p,:));
        title(sprintf('Tree %d leaf',p))
    end  
end

figure;plot_toydata(data_train);
hold on
plot(test_point(1,1), test_point(1,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.9 .5 .5], 'MarkerEdgeColor','k');
plot(test_point(2,1), test_point(2,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.5 .9 .5], 'MarkerEdgeColor','k');
plot(test_point(3,1), test_point(3,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.5 .5 .9], 'MarkerEdgeColor','k');
plot(test_point(4,1), test_point(4,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.9 .5 .5], 'MarkerEdgeColor','k');
title('Four Test Points Class Labels')

figure; 
for j = 1:4
    subplot(2,2,j);bar(p_rf_sum(j,:));
    title(sprintf('Final prediction of Point %d',j))
end



%% Q6 Test on the data_test
% Test on the dense 2D grid data, and visualise the results ... 
dense_leaves = testTrees_fastLinear(data_test,trees_linear);

% Get the probability of each data_test point from all TEN trees
for L = 1:length(data_test(:,1))
   p_rf_dense = trees(1).prob(dense_leaves(L,:),:); 
   p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end

% Visualise the 
fig = visualise(data_train,p_rf_dense_sum,[],0);
title('Visualise of Test Data Classification with Colour Encoded')






 


%% % Change the RF parameter values and evaluate ...
%Q7 Try different parameters of RF and see the effects
%% the number of trees
T_num = [1 5 10 20 50]; % give a set of tree numbers to test

init;
% Select dataset
[data_train, data_test] = getData('Toy_Spiral');

% Set the random forest parameters for instance, 
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied
for n = 1:length(T_num)   
    param.num = T_num(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(T_num)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test)
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense);
    end
    fig = visualise(data_train,p_rf_dense_sum,[],0);
    title(sprintf('Visualise of Test Data Classification of %d Trees',T_num(k)))
end

%% the depth of trees
init;
% Select dataset
[data_train, data_test] = getData('Toy_Spiral');

T_depth = [10 15]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.num = 10;        % trees num
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied
for n = 1:length(T_depth)   
    param.depth = T_depth(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(T_depth)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test)
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = sum(p_rf_dense)/length(diff_t{k});
    end
    fig = visualise(data_train,p_rf_dense_sum,[],0);
    title(sprintf('Visualise of Test Data Classification of depth %d',T_depth(k)))
end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% experiment with Caltech101 dataset for image categorisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init;

% Select dataset
[data_train, data_test] = getData('Caltech');

% we do bag-of-words technique to convert images to vectors (histogram of codewords)
% Set 'showImg' in getData.m to 0 to stop displaying training and testing images and their feature vectors
% close all;
%% Default settings
% Set the random forest parameters ...
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% Train Random Forest ...
figure;
trees = growTrees(data_train,param);
% Evaluate/Test Random Forest ...
dense_leaves = testTrees_fast(data_test,trees); 
% Get the probability of each data_test point from all TEN trees
for L = 1:length(data_test(:,1))
   p_rf_dense = trees(1).prob(dense_leaves(L,:),:); 
   p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end

% Visualise the 
% fig = visualise(data_train,p_rf_dense_sum,[],0);
% title('Visualise of Test Data Classification with Colour Encoded')
% show accuracy and confusion matrix ...
p_rf=p_rf_dense_sum;
figure;
confus_script

%% the number of trees
T_num = [1 5 10 20 50]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied

for n = 1:length(T_num)   
    figure;
    param.num = T_num(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(T_num)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
    figure;
    confus_script
end


%% the depth of trees 
% error when running depth of 20!!!!!!!!!!!
depth = [ 5 10 20 ]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.num = 10;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied

for n = 1:length(depth)   
    figure;
    param.depth = depth(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(depth)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
    figure;
    confus_script
end


%% Random forest codebook
[rf_data_train, rf_data_test]=rf_code('Caltech');
trees = growTrees(rf_data_train,param);
dense_leaves = testTrees_fast(rf_data_test,trees);
for L = 1:length(rf_data_test(:,1))
p_rf_dense = trees(1).prob(dense_leaves(L,:),:);
p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end
p_rf=p_rf_dense_sum;
figure;
confus_script


%%
clearvars -except data_test data_train Type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% random forest codebook for Caltech101 image categorisation
% .....



