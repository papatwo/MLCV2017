init;
[data_train, data_test] = getData('Toy_Spiral'); 

%% to create  bags of dataset 
% 4 bags without replacement
group = crossvalind('Kfold', length(data_train(:,1)), 4);
for i=1:4
    bag_re{i} = data_train(group==i,:); 
end

% bags with replacement of 63.2% unique data
frac = 1 - 1/exp(1); % Bootstrap sampling fraction: 1 - 1/e (63.2%)
[N,D]=size(data_train); 
idx = randsample(N,ceil(N*frac),1); 
%%for test and test git
