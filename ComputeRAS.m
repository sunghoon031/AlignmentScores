function [ RAS ] = ComputeRAS( R_input, R_gt )

    n_iterations = 10;
    
    n_total = size(R_gt, 3);
    R_samples = cell(1, n_total);
    for i = 1:n_total
        R_samples{i} = R_input(:,:,i)*R_gt(:,:,i)';
    end
    R_avg = RotationAveraging_TLUD(R_samples, n_iterations);
    
    errs = nan(1, n_total);
    for i = 1:n_total
        R_diff = (R_avg*R_gt(:,:,i))'*R_input(:,:,i);
        errs(i) = abs(acosd((trace(R_diff)-1)/2));
    end
    
    cumulative_relative_freq = [];
    thrs = 0.1:0.1:10;
    for thr = thrs
        cumulative_relative_freq(end+1) = sum(errs<thr)/(length(errs)*length(thrs));
    end
    RAS = sum(cumulative_relative_freq);

end

