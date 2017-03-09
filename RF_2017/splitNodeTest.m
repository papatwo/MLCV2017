function [node,nodeL,nodeR] = splitNodeTest(data,node,param,mode)
% Split node
% Param struct includes: 1.num; 2.depth; 3.splitNum; 4.split; 5.sp_func
% 6.threshold
% sp_func: axis, linear, 

visualise = 1;

% Initilise child nodes
iter = param.splitNum;
nodeL = struct('idx',[],'t',nan,'dim',0,'prob',[]);
nodeR = struct('idx',[],'t',nan,'dim',0,'prob',[]);

if length(node.idx) <= 5 % make this node a leaf if has less than 5 data points
    node.t = nan;
    node.dim = 0;
    return;
end

idx = node.idx;
data = data(idx,:);
[N,D] = size(data);
ig_best = -inf; % Initialise best information gain
idx_best = [];
   
for n = 1:iter
    
    % Split function - Modify here and try other types of split function
    if mode==1
        % axis-aligned split function:
        dim = randi(D-1); % Pick one random dimension -- pick x-axis or y-axis
        d_min = single(min(data(:,dim))) + eps; % Find the data range of this dimension
        d_max = single(max(data(:,dim))) - eps;
        
        if ~isfield(param,'threshold')
            t = d_min + rand*((d_max-d_min)); % Pick a random value within the range as threshold
        else
            t = param.threshold; % set a fixed threshold
        end
        idx_ = data(:,dim) < t;
        
        ig = getIG(data,idx_); % Calculate information gain: the point with max ig which satisfies the threshold
        % the info gain of all data at the left and right children nodes
    
        
    elseif mode == 2
        % linear split function
        if ~isfield(param,'threshold')
            t = 0; % Pick a random value within the range as threshold
        else
            t = param.threshold; % set a fixed threshold
        end
        r1= randi(D); %class number
        r2= randi(D); %class number
        w= randn(3, 1); %
        idx_ = [data(:, [r1 r2]), ones(N, 1)]*w < t;
        ig = getIG(data,idx_);
    end
    dim=3;
    [node, ig_best, idx_best] = updateIG(node,ig_best,ig,t,idx_,dim,idx_best);    
end

nodeL.idx = idx(idx_best);
nodeR.idx = idx(~idx_best);

end

function ig = getIG(data,idx) % Information Gain - the 'purity' of data labels in both child nodes after split. The higher the purer.
L = data(idx);
R = data(~idx);
H = getE(data);
HL = getE(L);
HR = getE(R);
ig = H - sum(idx)/length(idx)*HL - sum(~idx)/length(idx)*HR; 
end

function H = getE(X) % Entropy of all datapoints at the node!!!!
cdist= histc(X(:,1:end), unique(X(:,end))) + 1;
cdist= cdist/sum(cdist);
cdist= cdist .* log(cdist);
H = -sum(cdist);
end

function [node, ig_best, idx_best] = updateIG(node,ig_best,ig,t,idx,dim,idx_best) % Update information gain
if ig > ig_best
    ig_best = ig;
    if mode ==1
        node.t = t;
        node.dim = dim;
    elseif mode ==2
        node.dim=[r1,r2];
        node.t=w;
    end
    idx_best = idx;
else
    idx_best = idx_best;
end
end