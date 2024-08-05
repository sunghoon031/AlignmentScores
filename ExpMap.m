function out = ExpMap(in)
    angle = norm(in);
    if (angle == 0)
        out = eye(3);
        return;
    end
    axis = in/angle;
    ss = SkewSymmetricMatrix(axis);
    R = eye(3)+ss*sin(angle)+ss^2*(1-cos(angle));
    out = R;
end