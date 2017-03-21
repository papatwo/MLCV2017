param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';
% axis-aligned


% figure;
% test on different RF (num of trees)
% axis-aligned
for j=1:5
tree_A = growTrees(data_train,param);
dense_leaves = testTrees_fast(data_test,tree_A);
dense_leaves(dense_leaves==0)=1;
for L = 1:length(data_test(:,1))
    p_rf_dense =tree_A(1).prob(dense_leaves(L,:),:);
    p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end
p_rf=p_rf_dense_sum;
% subplot(1,2,1);
confus_script
accuracyA(j)=accuracy_rf;
end
% linear
for j=1:5
tree_L = growTreesLinear(data_train,param);
dense_leaves = testTrees_fastLinear(data_test,tree_L);
dense_leaves(dense_leaves==0)=1;
for L = 1:length(data_test(:,1))
    p_rf_dense =tree_L(1).prob(dense_leaves(L,:),:);
    p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end
p_rf=p_rf_dense_sum;
% subplot(1,2,2);
confus_script
accuracyL(j)=accuracy_rf;
end