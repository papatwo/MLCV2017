function d2 = euclid(x1,x2)
% Function to compute euclidean distance between
% two N-dimensional vectors, x1 and x2

N = max(size(x1));
if (N ~= max(size(x2)))
	disp('Error in euclid()');
	return;
end

d = x1 - x2;

d2 = sqrt(sum(d.^2));