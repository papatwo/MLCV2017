function [match, confidence, dist, ratio] = knn_match(featureA, featureB, threshold)
% threshold usually around 0.7-0.8
dist = zeros(size(featureA,1),size(featureB,1));
match = [];
confidence = [];
k = 1;
for i = 1:size(featureA,1)
    for j = 1:size(featureB,1)
        dist(i,j) = sqrt(sum((featureA(i,:)-featureB(j,:)).^2));
    end
    [dist_sort, idx_sort] = sort(dist(i,:), 'ascend'); 
    % Nearest Neighbor Distance Ratio Test
    ratio = dist_sort(1)/dist_sort(2);
    if ratio < threshold % find threshold by assessing max and min value of ratio
        match(k,1) = i; % match point for img A
        match(k,2) = idx_sort(1); % match point for img B
        confidence(k) = ratio; 
        k = k+1;
    end
end
	
% a=intere_ptA(:,6)'
% figure(5);hold on;plot(a(1,1),a(1,2),'rs')
% b=intere_ptB(:,19)'
% figure(7);hold on;plot(b(1),b(2),'rs')
end