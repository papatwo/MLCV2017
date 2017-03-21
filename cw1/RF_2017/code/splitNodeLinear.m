function [node,nodeL,nodeR] = splitNodeLinear(data,node,param)
% Split node
% Param struct includes: 1.num; 2.depth; 3.splitNum; 4.split; 5.sp_func
% 6.threshold
% sp_func: axis, linear, 


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

        % linear split function

        r1= randi(D-1); %class number
        r2= randi(D-1); %class number
        w= randn(3, 1); %
        t=w;
        dim=[r1,r2];
        idx_ = [data(:, [r1 r2]), ones(N, 1)]*w < 0;
        ig = getIG(data,idx_);

    
    [node, ig_best, idx_best] = updateIG(node,ig_best,ig,t,idx_,dim,idx_best);    
end

nodeL.idx = idx(idx_best);
nodeR.idx = idx(~idx_best);

end

function ig = getIG(data,idx) % Information Gain - the 'purity' of data labels in both child nodes after split. The higher the purer.
data=data(:,end);
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
        node.t = t;
        node.dim = dim;
    idx_best = idx;
else
    idx_best = idx_best;
end
end