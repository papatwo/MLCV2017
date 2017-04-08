function H = get_homography(corre1,corre2)
% input points vector must be 2xN size and at least contain 4 points

% x1 = corre1(1,:);
% y1 = corre1(2,:);
% x2 = corre2(1,:);
% y2 = corre2(2,:);
n = size(corre1,2);
A = zeros(2*n,9);

for i = 1:n % for loop to create A matrix where Ah=0
    A(2*i-1,:) = [corre1(:,i)',1,0,0,0,-corre1(:,i)'.*corre2(1,i),-corre2(1,i) ];
    A(2*i,:) = [0,0,0,corre1(:,i)',1,-corre1(:,i)'.*corre2(2,i),-corre2(2,i) ];
end

if n>=4 % there are right no. equations to solve 8 unknows in H
    [~,~,V] = svd(A,0);
    H = reshape(V(:,end),[3 3])';
else
    error('Not enough points');
end
H = H./H(end,end);
end

