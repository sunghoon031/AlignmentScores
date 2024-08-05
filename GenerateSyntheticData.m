function [xyz_gt, R_gt, xyz_input, R_input] = GenerateSyntheticData(n_total, n_outlier, sigma_xyz, sigma_R)

    % Ground truth points are inside 1x1x1 cube centered at the origin.
    xyz_gt = nan(3, n_total);
    R_gt = nan(3,3, n_total);

    for i = 1:n_total
        xyz_gt(:,i) = rand(3,1)-0.5;
        R_gt(:,:,i) = RandomRotationMatrix;
    end


    % Outlier points are inside 10x10x10 cube centered at the origin.
    xyz_outliers = nan(3, n_total);
    R_outliers = nan(3,3, n_total);

    for i = 1:n_total
        xyz_outliers(:,i) = (rand(3,1)-0.5)*10;
        R_outliers(:,:,i) = RandomRotationMatrix;
    end

    xyz_input = [xyz_gt(:,1:n_total-n_outlier), xyz_outliers(:, 1:n_outlier)];
    xyz_input = xyz_input + normrnd(0, sigma_xyz, size(xyz_input));

    R_input = nan(size(R_gt)); 
    R_transform = RandomRotationMatrix; % GT world to EST world
    s_transform = rand(1)*10;
    t_transform = rand(3,1)*100;

    for i = 1:n_total
        R_input(:,:,i) = RandomRotation(abs(normrnd(0, sigma_R)))*R_gt(:,:,i);
        R_input(:,:,i) = R_transform*R_input(:,:,i);

        xyz_input(:,i) = s_transform*R_transform*xyz_input(:,i) + t_transform;
    end
    R_input(:,:,n_total-n_outlier+1:end) = R_outliers(:,:,1:n_outlier);



end