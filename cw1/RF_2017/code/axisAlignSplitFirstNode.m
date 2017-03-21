% splitting the first node
bestig= -100;
iter = param.splitNum;
data_train = bag{1,1}; % use one subbag to split one tree
for n = 1:iter
    visualise = 1;
    % axis-aligned split function:
    dim = randi(D-1); % Pick one random dimension -- pick x-axis or y-axis
    d_min = single(min(data_train(:,dim))) + eps; % Find the data range of this dimension
    d_max = single(max(data_train(:,dim))) - eps;
    t = d_min + rand*((d_max-d_min)); % Pick a random value within the range as threshold
    
    idx_ = data_train(:,dim) < t; % sort out data points go LEFT by the threshold
    l_data = data_train(idx_,:); % data go LEFT branch
    r_data = data_train(~idx_,:); % data go RIGHT branch
    
    ig = getIG(data_train, idx_);
    
%     % calculate entropy   
%     H = getE(data_train); % before split
%     HL = getE(l_data); %after split - L
%     HR = getE(r_data); % - R
%     % calculate information gain and update the best ig
%     ig = H - (sum(idx_)/length(idx_)*HL + sum(~idx)/length(idx)*HR);
    
    if bestig < ig
        bestig = ig;
        t_best = t;
        dim_best = dim;
        idx_best = idx_;
    end
    
    if visualise
        visualise_splitfuncTest(idx_,data_train,dim,t,ig,0);
        pause();
    end  
end
disp('showing best')
% figure; % visualise the best split
visualise_splitfuncTest(idx_best,data_train,dim_best,t_best,bestig,0);