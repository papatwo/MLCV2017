splitNum = [3 5 8 10 15]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.num = 10;        % trees num
param.depth = 5;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% Axis!!!!!!!!!!!
% grow corresponding RF with respect to the num of trees varied
for n = 1:length(splitNum)   
    param.num = splitNum(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(splitNum)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    dense_leaves(dense_leaves==0)=1;
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
    subplot(2,3,k);
    confus_script
end
