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
