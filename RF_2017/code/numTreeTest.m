T_num = [5 10 20 50 100 200]; % give a set of tree numbers to test

% init;
% % Select dataset
% [data_train, data_test] = getData('Toy_Spiral');

% Axis!!!!!!!!!!!
% grow corresponding RF with respect to the num of trees varied
for n = 1:length(T_num)   
    param.num = T_num(n);         % Number of trees
    diff_tA{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(T_num)
    subplot(2,3,k);
    trees_axis = diff_tA{k};
    testRF_Axis;
    str=sprintf('%d Trees RF',T_num(k));
    title(str)
end

% Linear!!!!!!!!!!!
% grow corresponding RF with respect to the num of trees varied
for n = 1:length(T_num)   
    param.num = T_num(n);         % Number of trees
    diff_tL{n} = growTreesLinear(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(T_num)
    subplot(2,3,k);
    trees_linear = diff_tL{k};
    testRF_Linear;
    str=sprintf('%d Trees RF',T_num(k));
    title(str)
end