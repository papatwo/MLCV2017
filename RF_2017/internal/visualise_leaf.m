function visualise_leaf(trees)
for L = 1:9
    try
        subplot(3,3,L);
        bar(trees(1).leaf(L).prob);
        title(sprintf('Leaf node %d',L))
        axis([0.5 3.5 0 1]);
    end
end
end