function [errors, s_final, R_final, t_final] = ModifiedPCR99(xyz_gt,xyz_est, thr1, n_hypo)
    
    % Modified from https://github.com/sunghoon031/PCR-99/blob/main/PCR99b.m
 
    n = size(xyz_gt, 2);
    
    m = max(4, round(n*0.1));
    
   
    
    log_ratio_mat = nan(n,n);
    
    min_mth_error = inf;
    c = 0;
    
    for it = 1:10^9

        ijk = randperm(n, 3);
        i = ijk(1);
        j = ijk(2);
        k = ijk(3);


        if (isnan(log_ratio_mat(i,j)))
            d_gt_ij = norm(xyz_gt(:,i) - xyz_gt(:,j));
            d_est_ij = norm(xyz_est(:,i) - xyz_est(:,j));

            log_ratio_mat(i,j) = log(d_est_ij/d_gt_ij);
            log_ratio_mat(j,i) = log_ratio_mat(i,j);
        end
        if (isnan(log_ratio_mat(j,k)))
            d_gt_jk = norm(xyz_gt(:,j) - xyz_gt(:,k));
            d_est_jk = norm(xyz_est(:,j) - xyz_est(:,k));

            log_ratio_mat(j,k) = log(d_est_jk/d_gt_jk);
            log_ratio_mat(k,j) = log_ratio_mat(j,k);
        end
        if (isnan(log_ratio_mat(k,i)))
            d_gt_ki = norm(xyz_gt(:,k) - xyz_gt(:,i));
            d_est_ki = norm(xyz_est(:,k) - xyz_est(:,i));

            log_ratio_mat(k,i) = log(d_est_ki/d_gt_ki);
            log_ratio_mat(i,k) = log_ratio_mat(k,i);
        end

        log_ratio_ij = log_ratio_mat(i,j);
        log_ratio_jk = log_ratio_mat(j,k);
        log_ratio_ki = log_ratio_mat(k,i);

        e1 = abs(log_ratio_ij - log_ratio_jk);
        e2 = abs(log_ratio_jk - log_ratio_ki);
        e3 = abs(log_ratio_ki - log_ratio_ij);

        if (e1 > thr1 || e2 > thr1 || e3 > thr1)
            continue;
        end

        c = c + 1;

        A = xyz_gt(:,[i, j, k]);
        B = xyz_est(:,[i, j, k]);


        [scale, R, t] = sRt_from_3points(A,B);
        
        
        errs = sqrt(sum(((xyz_est-t)/scale - R*xyz_gt).^2,1));
        
        errs_sorted = sort(errs);
        
        mth_error = errs_sorted(m);
        
        if (mth_error < min_mth_error)
            min_mth_error = mth_error;
            R_final = R;
            t_final = t;
            s_final = scale;
        end
        
        if (c==n_hypo)
            errors = sqrt(sum(((xyz_est-t_final)/s_final - R_final*xyz_gt).^2,1));
            break;
        end
    end

end