splitNum = [3 5 8 10 15]; % give a set of tree numbers to test
% Set the random forest parameters for instance, 
param.num = 10;        % trees num
param.depth = 5;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% Axis!!!!!!!!!!!
% grow corresponding RF with respect to the num of trees varied
for n = 1:length(splitNum)   
    param.num = splitNum(n);         % Number of trees
    diff_tA{n} = growTrees(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(splitNum)    
    subplot(2,3,k);
    trees_axis = diff_tA{k};
    testRF_Axis;
    str=sprintf('%d Split RF',splitNum(k));
    title(str)
end

% Linear!!!!!!!!!!!
% grow corresponding RF with respect to the num of trees varied
for n = 1:length(splitNum)   
    param.num = splitNum(n);         % Number of trees
    diff_tL{n} = growTreesLinear(data_train,param);
end
figure;
% test on different RF (num of trees)
for k = 1:length(splitNum)   
    subplot(2,3,k);
    trees_linear = diff_tL{k};
    testRF_Linear;
    str=sprintf('%d Split RF',splitNum(k));
    title(str)
end