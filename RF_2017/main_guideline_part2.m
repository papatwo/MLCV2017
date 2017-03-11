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
subplot(2,2,1)
kcode_RF_default;


%% the number of trees
kcode_RF_numTree;
accT=mean(accuracy,2)';
numT=[5 10 20 50 100 200 300]; % for plotting
% plot(numT,acc,'*-r')
%% the depth of trees 
kcode_RF_depth;
accD=  [0.2413 0.4133 0.4867 0.4720 0.4493]; % tested on uni PC as laptop nor powerful enough
numD = [2 5 8 10 15];
%% number of split
kcode_RF_splitNum;
accS=mean(accuracy,2)';
numS=[3 5 8 10 15 20 ]; % for plotting

subplot(2,2,2)
plot(numT,accT,'*-b');
subplot(2,2,3)
plot(numD,accD,'*-b');
subplot(2,2,4);
plot(numS,accS,'*-b');
%% Split function type
kcode_RF_wl;
accA=mean(accuracyA,2)';
accL=mean(accuracyL,2)';
figure;subplot(2,2,1)
plot([1:5],accuracyA,'*-r');
hold on
plot([ 1:5],accuracyL,'*-b');

%% Best kmeans codebook tree
param.num = 200;         % Number of trees
param.depth = 8;        % trees depth
param.splitNum = 50;     % Number of split functions to try
param.split = 'IG';
% Train Random Forest ...
figure;
 for i=1:5
tic
trees = growTrees(data_train,param);
train_t=toc;
% Evaluate/Test Random Forest ...
% dense_leaves = testTrees_fastLinear([data_train(:,1:end-1),zeros(size(data_test,1),1)],trees); 

dense_leaves = testTrees_fast([data_test(:,1:end-1),zeros(size(data_test,1),1)],trees); 
dense_leaves(dense_leaves==0)=1;
% Get the probability of each data_test point from all TEN trees
tic
for L = 1:length(data_test(:,1))
   p_rf_dense = trees(1).prob(dense_leaves(L,:),:); 
   p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end
test_t=toc;
% Visualise the 
% fig = visualise(data_train,p_rf_dense_sum,[],0);
% title('Visualise of Test Data Classification with Colour Encoded')
% show accuracy and confusion matrix ...
p_rf=p_rf_dense_sum;
subplot(2,3,i)
confus_script
 end

%% Random forest codebook

[rf_data_train, rf_data_test,booktime,codebook]=rf_code('Caltech');

param.num = 200;         % Number of trees
param.depth = 8;        % trees depth
param.splitNum = 50;     % Number of split functions to try
param.split = 'IG';
trees = growTrees(rf_data_train,param);

dense_leaves = testTrees_fast(rf_data_test(:,1:end-1),trees);
dense_leaves(dense_leaves==0)=1;
for L = 1:length(rf_data_test(:,1))
p_rf_dense = trees(1).prob(dense_leaves(L,:),:);
p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end
p_rf=p_rf_dense_sum;
figure;
data_test=rf_data_test;
confus_script

% tried different vocab size used for RF codebook and testing by the
% best settings obtained from q32 to see the effects on accuracy
figure;
subplot(2,2,1);plot([32 76 125 231],[0.613 0.667 0.68 0.653],'*-b')
subplot(2,2,2);plot([3 5 10 20 50],[0.66 0.727 0.687 0.633 0.58],'*-b')
subplot(2,2,3);plot([3 5 10 15],[0.52 0.687 0.687 0.653],'*-b')
%% RF CLASSIFIER param tuning
[rf_data_train, rf_data_test,booktime,codebook]=rf_code('Caltech'); %currently is the best RF param for codebook
data_train=rf_data_train;
data_test=rf_data_test;
kcode_RF_default;
title('RF codebook Classifier Test')
%% the number of trees
kcode_RF_numTree;
accT=mean(accuracy,2)';
numT=[5 10 20 50 100 200 300];
%% the depth of trees 
kcode_RF_depth;
accD=  [0.2467 0.4267 0.4067 0.46 0.4467]; % tested on uni PC as laptop nor powerful enough
numD = [2 5 8 10 15];
%% number of split
kcode_RF_splitNum;
accS=mean(accuracy,2)';
numS=[3 5 8 10 15 20 ]; % for plotting

subplot(2,2,2)
plot(numT,accT,'*-b');
subplot(2,2,3)
plot(numD,accD,'*-b');
subplot(2,2,4);
plot(numS,accS,'*-b');
%%
clearvars -except data_test data_train Type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% random forest codebook for Caltech101 image categorisation
% .....



