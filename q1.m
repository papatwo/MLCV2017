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
    bag{i} = data_train(randsample(N,ceil(N*frac),1),:); % 4 bags of samples with replacement 
    figure; plot_toydata(bag{i}) 
end




%% grow trees 

%% Set the random forest parameters for instance, 
param.num = 10;         % Number of trees
param.depth = 5;        % trees depth
param.splitNum = 3;     % Number of split functions to try (what's this exactly mean???
param.split = 'IG';     % Currently support 'information gain' only

% grow all trees
trees = growTrees(data_train,param);
% prob field in trees struct: only the prob of first tree was calculated.
% to get all ten trees prob go to growTrees.m to line 47 and 49 to modify
% index of trees.....

%% Use Karpathy RF library
opts= struct;
opts.depth= 5;
opts.numTrees= 10;
opts.numSplits= 3;
opts.verbose= false;
opts.classifierID= 2; % weak learners to use. Can be an array for mix of weak learners too

X= bsxfun(@rdivide, bsxfun(@minus, data_train(:,1:2), mean(data_train(:,1:2))), var(data_train(:,1:2)));
Y=data_train(:,3);
% classify
tic;
m= forestTrain(X, Y, opts);
timetrain= toc;
tic;
yhatTrain = forestTest(m, X);
timetest= toc;
bar(hist(Y,unique(Y)));
bar(hist(YR,unique(YR)));
bar(hist(YL,unique(YL)));
%%%change input to subbags data


%%

% Look at classifier distribution for fun, to see what classifiers were
% chosen at split nodes and how often
fprintf('Classifier distributions:\n');
classifierDist= zeros(1, 4);
unused= 0;
for i=1:length(m.treeModels)
    for j=1:length(m.treeModels{i}.weakModels)
        cc= m.treeModels{i}.weakModels{j}.classifierID;
        if cc>1 %otherwise no classifier was used at that node
            classifierDist(cc)= classifierDist(cc) + 1;
        else
            unused= unused+1;
        end
    end
end
fprintf('%d nodes were empty and had no classifier.\n', unused);
for i=1:4
    fprintf('Classifier with id=%d was used at %d nodes.\n', i, classifierDist(i));
end

%% plot results
xrange = [-1.5 1.5];
yrange = [-1.5 1.5];
inc = 0.02;
[x, y] = meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2));
image_size = size(x);
xy = [x(:) y(:)];

[yhat, ysoft] = forestTest(m, xy);
decmap= reshape(ysoft, [image_size 3]);
decmaphard= reshape(yhat, image_size);

subplot(121);
imagesc(xrange,yrange,decmaphard);
hold on;
set(gca,'ydir','normal');
cmap = [1 0.8 0.8; 0.95 1 0.95; 0.9 0.9 1];
colormap(cmap);
plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
hold off;
title(sprintf('%d trees, Train time: %.2fs, Test time: %.2fs\n', opts.numTrees, timetrain, timetest));

subplot(122);
imagesc(xrange,yrange,decmap);
hold on;
set(gca,'ydir','normal');
plot(X(Y==1,1), X(Y==1,2), 'o', 'MarkerFaceColor', [.9 .3 .3], 'MarkerEdgeColor','k');
plot(X(Y==2,1), X(Y==2,2), 'o', 'MarkerFaceColor', [.3 .9 .3], 'MarkerEdgeColor','k');
plot(X(Y==3,1), X(Y==3,2), 'o', 'MarkerFaceColor', [.3 .3 .9], 'MarkerEdgeColor','k');
hold off;