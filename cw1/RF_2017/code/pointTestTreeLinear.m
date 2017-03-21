test_point = [-.5 -.7; .4 .3; -.7 .4; .5 -.5];
figure;plot_toydata(data_train);
plot(test_point(:,1), test_point(:,2), 's', 'MarkerSize',20, 'MarkerFaceColor', [.9 .9 .9], 'MarkerEdgeColor','k');
title('Visualisation of Test Points');

for n=1:length(test_point)%4 % how to get all 22801 test point prob
    leaves = testTreesLinear([test_point(n,:) 0],trees_linear);
    [zr,~] = find(leaves==0);
    leaves(leaves==0) = 1;
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
