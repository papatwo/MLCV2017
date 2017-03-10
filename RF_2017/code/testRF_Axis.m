% Test on the dense 2D grid data, and visualise the results ... 
dense_leaves = testTrees_fast(data_test,trees_axis);
% Get the probability of each data_test point from all TEN trees
for L = 1:length(data_test(:,1))
   p_rf_dense = trees_axis(1).prob(dense_leaves(L,:),:); 
   p_rf_dense_sum(L,:) = mean(p_rf_dense,1);
end

% % Visualise the 
% fig = visualise(data_train,p_rf_dense_sum,[],0);
% title('Visualise of Test Data Classification with Colour Encoded')

c1=[];c2=[];c3=[];
for i=1:length(p_rf_dense_sum)
     [a b]=max(p_rf_dense_sum(i,:));
     if b==1
         c1=[c1 i];
     elseif b==2
         c2=[c2 i];
     else
         c3=[c3 i];
     end
end
plot_toydata(data_train);
scatter(data_test(c1,1),data_test(c1,2),'.r');
scatter(data_test(c2,1),data_test(c2,2),'.g');
scatter(data_test(c3,1),data_test(c3,2),'.b');
plot_toydata(data_train);