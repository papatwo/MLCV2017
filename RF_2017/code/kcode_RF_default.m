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
% figure;
confus_script