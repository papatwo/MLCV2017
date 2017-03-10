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
splitNum=[1 3 5 8 10 ];
for i=1:length(splitNum)
    param.splitNum=splitNum(i);
    axisAlignSplitFirstNode;
    disp('finish')
end


%% Q3 Grow ONE complete tree
T=1;
param.num=1;
tree_axis=growTrees(data_train,param);
trees = tree_axis;
% Store class distribution in the leaf nodes.
makeLeaf;
% Visualise the leaf node class distribution
figure;
visualise_leaf(tree_axis)
% every time train a tree will have different number of leaf node

%% linear grow one tree
T = 1;
param.num=1;
tree_linear=growTreesLinear(data_train,param);
% Store class distribution in the leaf nodes.
makeLeaf;
% Visualise the leaf node class distribution
figure;
visualise_leaf

%% Q4 Grow all 10 trees
% axis-aligned
param.num=10;
trees_axis = growTrees(data_train,param); 
% linear
param.num=10;
trees_linear=growTreesLinear(data_train,param);

%% Q5 Testing data points
% Evaluate/Test Random Forest grab the few data points and evaluate them
% one by one by the leant RF

% axis-aligned
pointTestTreeAxis;
% linear
pointTestTreeLinear;

%% Q6 Test on the data_test

% axis-aligned
testRF_Axis;

% linear
testRF_Linear;

%% % Change the RF parameter values and evaluate ...
%Q7 Try different parameters of RF and see the effects

% Set the random forest parameters by default 
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only
% the number of trees
numTreeTest;

% the depth of trees
depthTreeTest;

% the num of split
numSplitTest;

% grow one optimal tree for each method
param.num = 200;         % Number of trees
param.depth = 10;        % trees depth
param.splitNum = 8;     % Number of split functions to try
param.split = 'IG';

b_treeA = growTrees(data_train,param);
b_treeL = growTreesLinear(data_train,param);

figure;
trees_axis = b_treeA;
testRF_Axis;

figure;
trees_linear = b_treeL;
testRF_Linear;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% experiment with Caltech101 dataset for image categorisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init;

% Select dataset
[data_train, data_test,booktime] = getData('Caltech');

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
trees = growTrees(data_train,param);
% Evaluate/Test Random Forest ...
dense_leaves = testTrees_fast([data_test(:,1:end-1),zeros(size(data_test,1),1)],trees); 
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
T_num = [5 10 20 50 100 200]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied

for n = 1:length(T_num)   
    param.num = T_num(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(T_num)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
    subplot(2,3,k);
    confus_script
%     str=sprintf('%d Trees RF',T_num(k));
%     title(str)
end


%% the depth of trees 
% error when running depth of 20!!!!!!!!!!!
depth = [ 2 5 8 10 15 50 ]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.num = 10;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied

for n = 1:length(depth)   
    param.depth = depth(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(depth)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
    subplot(2,3,k);
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



