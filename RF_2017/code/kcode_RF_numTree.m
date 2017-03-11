T_num = [5 10 20 50 100 200 300]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow corresponding RF with respect to the num of trees varied
for j = 1:5
for n = 1:length(T_num)   
    param.num = T_num(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end
figure;
% test on different RF (num of trees)

for k = 1:length(T_num)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    dense_leaves(dense_leaves==0)=1;
    for L = 1:length(data_test(:,1))
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
    end
    p_rf=p_rf_dense_sum;
%     subplot(2,3,k);
    confus_script
    accuracy (k,j) = accuracy_rf; 
%     str=sprintf('%d Trees RF',T_num(k));
%     title(str)
end
end
