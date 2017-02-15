init;
[data_train, data_test] = getData('Toy_Spiral'); 

%% to create  bags of dataset 
% 4 bags without replacement
% group = crossvalind('Kfold', length(data_train(:,1)), 4);
% for i=1:4
%     bag_uni{i} = data_train(group==i,:); 
% end

% bags with replacement
frac = 1 - 1/exp(1); % Bootstrap sampling fraction: 1 - 1/e (63.2%)
[N,D]=size(data_train);
bag={};
for i=1:4
    bag{i} = data_train(randsample(N,ceil(N*frac),1)); % 4 bags of samples with replacement  
end

%% grow trees
% Set the random forest parameters for instance, 
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try
param.split = 'IG';     % Currently support 'information gain' only

% grow all trees
trees = growTrees(data_train,param);
