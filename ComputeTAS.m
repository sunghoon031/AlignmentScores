function [TAS] = ComputeTAS(xyz_gt,xyz_input)

    thr1 = 0.1;
    n_hypo = 1000;
    [errs, ~,~,~] = ModifiedPCR99(xyz_gt,xyz_input, thr1, n_hypo);

    thr = GetThreshold(xyz_gt);

    cumulative_relative_freq = [];
    thrs = (0.01:0.01:1)*thr;
    for thr = thrs
        cumulative_relative_freq(end+1) = sum(errs<thr)/(length(errs)*length(thrs));
    end
    TAS = sum(cumulative_relative_freq);

end

