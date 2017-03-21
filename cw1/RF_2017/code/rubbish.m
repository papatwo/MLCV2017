%% CHANGE SPLIT FUNCTION??? -----NO USE 


% splitting the first node
bestig= -100;
iter = param.splitNum;
data = bag{1,1}; % use one subbag to split one tree
[N,D]=size(data);
for n = 1:iter
    visualise = 0;
    % axis-aligned split function:
    r1 = randi(D-1);
    r2 = randi(D-1); % Pick one random dimension -- pick x-axis or y-axis
    w= randn(3, 1); %
    idx_ = [data(:, [r1 r2]), ones(N, 1)]*w < 0; 
    dim = [r1 r2];
    t = w;
    l_data = data_train(idx_,:); % data go LEFT branch
    r_data = data_train(~idx_,:); % data go RIGHT branch
    
    ig = getIG(data_train, idx_);
    
    if bestig < ig
        bestig = ig;
        t_best = t;
        dim_best = dim;
        idx_best = idx_;
    end
    
%     if visualise
%         visualise_splitfunc(idx_,data_train,dim,t,bestig,0);
%         pause();
%     end    
end

figure; % visualise the best split
visualise_splitfuncTest(idx_best,data_train,dim_best,t_best,bestig,0);



% idx=1:length(bag{1}); 
% trees(T).node(1) = struct('idx',idx,'t',nan,'dim',-1,'prob',[]);
% 
% % Split Nodes
% for n = 1:2^(param.depth-1)-1
%     [trees(T).node(n),trees(T).node(n*2),trees(T).node(n*2+1)] = splitNodeLinear(data_train,trees(T).node(n),param);
% end




% grow corresponding RF with respect to the num of trees varied
for n = 1:length(T_depth)   
    param.depth = T_depth(n);         % Number of trees
    diff_t{n} = growTrees(data_train,param);
end

% test on different RF (num of trees)
for k = 1:length(T_depth)
    dense_leaves = testTrees_fast(data_test,diff_t{k});
    for L = 1:length(data_test)
        p_rf_dense = diff_t{k}(1).prob(dense_leaves(L,:),:);
        p_rf_dense_sum(L,:) = sum(p_rf_dense)/length(diff_t{k});
    end
    fig = visualise(data_train,p_rf_dense_sum,[],0);
    title(sprintf('Visualise of Test Data Classification of depth %d',T_depth(k)))
end
