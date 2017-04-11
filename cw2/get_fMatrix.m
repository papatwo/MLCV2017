function [F,e1,e2] = get_fMatrix(pt_l, pt_r)
% compute fundatmental matrix 
% (x_l)'F(x_r) = 0 where F is the 3x3 fundamental matrix
% F = [a b c; d e f; g h m]; 
% 9 components; Rank 2; 7 degrees of freedom
% F captures the relationship between the corresponding points in two views
% require at least 8 points (altho 9 unknows but its homogeneous system)

% input arg: 
% pt_l & pt_r are 2xN matching points matrix
% N at least be 8

% output:
%          F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
%          e1     - The epipole in image 1 such that F*e1 = 0
%          e2     - The epipole in image 2 such that F'*e2 = 0

if size(pt_l,2)<3
    pt_l = [pt_l ones(size(pt_l,1),1)]';
    pt_r = [pt_r ones(size(pt_r,1),1)]';
end
% Normalise each set of points so that the origin
% is at centroid and mean distance from origin is sqrt(2).
% normalise2dpts also ensures the scale parameter is 1.
[pt_l, T1] = normalise2dpts(pt_l);
[pt_r, T2] = normalise2dpts(pt_r);
% by expanding (pt_l)'*F*(pt_r), rewrite to A*f
A = [pt_r(1,:)'.*pt_l(1,:)'   pt_r(1,:)'.*pt_l(2,:)'  pt_r(1,:)' ...
    pt_r(2,:)'.*pt_l(1,:)'   pt_r(2,:)'.*pt_l(2,:)'  pt_r(2,:)' ...
    pt_l(1,:)'            pt_l(2,:)'            ones(size(pt_l,2),1) ]; 

[~,~,V] = svd(A,0);
F = reshape(V(:,end),[3,3])';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';
% Denormalise
F = T2'*F*T1;

[U,~,V] = svd(F,0);
e1 = hnormalise(V(:,3));
e2 = hnormalise(U(:,3));
% http://ece631web.groups.et.byu.net/Lectures/ECEn631%2013%20-%208%20Point%20Algorithm.pdf
end
