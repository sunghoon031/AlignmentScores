function R_out = RotationAveraging_TLUD(R_samples, nSteps)

    n_samples = length(R_samples);
    
    inliers = 1:n_samples;
    
    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_samples{i}(:);
    end                 

    
    thr = 0.5;
    min_cost = inf;
    for i = 1:n_samples
        vec_i = vectors_total(:,i);
        
        vec_diff = vectors_total-vec_i;
        vec_diff = vec_diff.^2;
        es = sqrt(sum(vec_diff,1));
        
        
        es(es>thr) = thr;
        cost = sum(es);
        
        if (cost < min_cost)
            min_cost = cost;
            inliers = find(es<thr);
        end
    
    end

    
    R_sum = zeros(3,3);
    for i = inliers
        R_sum = R_sum + R_samples{i};
    end
    [U,~,V] = svd(R_sum);
    R_out = U*V.';
    if (det(R_out) < 0)
        V(:,3) = -V(:,3);
        R_out = U*V.';
    end
    
    n_inliers = length(inliers);

    q_out = R2q(R_out);
    
    q_all = zeros(4,n_inliers);
    for i = 1:n_inliers
        q_all(:,i) = R2q(R_samples{inliers(i)});
    end


    q_all_inv_q_out = zeros(4, n_inliers);

    for j = 1:nSteps

        q_all_inv_q_out(1,:) = ... %scalar term of (q_all)*inv(q_out)
            -q_all(1,:).*q_out(1,:) - sum(q_all(2:4,:).*q_out(2:4,:),1); 

        q_all_inv_q_out(2:4,:) = ... % vector term of (q_all)*inv(q_out)
           +q_all(1,:).*q_out(2:4,:) - q_out(1,:).*q_all(2:4,:) ...
           + [q_all(3,:).*q_out(4,:) - q_all(4,:).*q_out(3,:);...
              q_all(4,:).*q_out(2,:) - q_all(2,:).*q_out(4,:);...
              q_all(2,:).*q_out(3,:) - q_all(3,:).*q_out(2,:)];


        sine_half = sqrt(sum(q_all_inv_q_out(2:4,:).^2, 1));
        theta = 2*atan2(sine_half, q_all_inv_q_out(1,:));
        theta(theta < -pi) = theta(theta < -pi) + 2*pi;
        theta(theta > pi) = theta(theta > pi) - 2*pi;

        q_all_inv_q_out(2:4,:) = q_all_inv_q_out(2:4,:).*sign(theta);
        theta = abs(theta);

        unit_v = q_all_inv_q_out(2:4,:)./sine_half;

        delta = sum(unit_v, 2)/sum(1./theta);
        delta_angle = norm(delta);

        unit_delta = delta/delta_angle;

        q_delta = zeros(4,1);
        q_delta(1) = cos(delta_angle/2);
        q_delta(2:4) = unit_delta*sin(delta_angle/2);

        q_out_ = q_out;
        q_out_(1,:) = ... %scalar term of (q_delta)*q_geo1
            q_delta(1).*q_out(1) - sum(q_delta(2:4).*q_out(2:4),1); 

        q_out_(2:4,:) = ... % vector term of (q_delta)*q_geo1
           q_delta(1).*q_out(2:4) + q_out(1,:).*q_delta(2:4) ...
           + [q_delta(3).*q_out(4) - q_delta(4,:).*q_out(3);...
              q_delta(4).*q_out(2) - q_delta(2,:).*q_out(4);...
              q_delta(2).*q_out(3) - q_delta(3,:).*q_out(2)];

        q_out = q_out_;  

        if (delta_angle < 0.001)
            break;
        end
    end

    R_out = q2R(q_out);

end
