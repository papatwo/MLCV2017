function HN = get_homography_norm(corre1,corre2)
% input points vector must be 2xN size and at least contain 4 points

x1 = corre1(1,:);
y1 = corre1(2,:);
x2 = corre2(1,:);
y2 = corre2(2,:);

n = size(corre1,2);
A = zeros(2*n,9);

S1 = sqrt(2)*n /sum(sqrt(power(x1-mean(x1),2) + power(y1-mean(y1),2)));

Norm1 = S1*[1, 0, -mean(x1);
    0, 1, -mean(y1);
    0, 0, 1/S1];

S2 = sqrt(2)*n /sum(sqrt(power(x2-mean(x2),2) + power(y2-mean(y2),2)));

Norm2 = S2*[1, 0, -mean(x2);
    0, 1, -mean(y2);
    0, 0, 1/S2];

corre1 = Norm1*[corre1; ones(1,n)];
corre2 = Norm2*[corre2; ones(1,n)];

corre1(3,:) = [];
corre2(3,:) = [];

for i = 1:n % for loop to create A matrix where Ah=0
    A(2*i-1,:) = [corre1(:,i)',1,0,0,0,-corre1(:,i)'.*corre2(1,i),-corre2(1,i) ];
    A(2*i,:) = [0,0,0,corre1(:,i)',1,-corre1(:,i)'.*corre2(2,i),-corre2(2,i) ];
end

if n>=4 % there are right no. equations to solve 8 unknows in H
    [~,~,V] = svd(A,0);
    HN = reshape(V(:,end),[3 3])';
else
    error('Not enough points');
end
HN = HN./HN(end,end);
HN = Norm2 / HN * Norm1;
end

