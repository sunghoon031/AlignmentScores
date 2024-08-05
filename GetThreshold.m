function thr = GetThreshold(xyz_gt)

    % Find the upper quartile of the closese-pair-distances.
    
    n = size(xyz_gt, 2);
    
    min_dists = zeros(1,n);
    
    for i = 1:n
        min_dist = inf;
        for j = 1:n
            if (i==j)
                continue;
            end
            
            d = norm(xyz_gt(:,i) - xyz_gt(:,j));
            if (d < min_dist)
                min_dist = d;
            end
        end
        
        min_dists(i) = min_dist;
    end
        
    dists = sort(min_dists);
    
    m = ceil(0.75*n);
    thr = dists(m);
        
end
