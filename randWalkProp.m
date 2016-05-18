function cand = randWalkProp(ro,breakpoints)
    % breakpoints contains t_i-1, t_i, t_i+1
    R = ro*(breakpoints(3)-breakpoints(1));
    eps = -R + (2*R)*rand;
    cand = breakpoints(2) + eps;    


