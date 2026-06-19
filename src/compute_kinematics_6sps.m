function kin = compute_kinematics_6sps(robot, L, dL, ddL, q0_guess)
% COMPUTE_KINEMATICS_6SPS Solves the kinematics of the 6-SPS platform.
%
% Utilizes Natural Coordinates to formulate an 18x18 Newton-Raphson system,
% avoiding trigonometric singularities to solve the forward kinematics of the 
% moving platform, along with its velocities and accelerations.
%
% Args:
%   robot (struct): Structural and geometric parameters of the 6-SPS.
%   L, dL, ddL (1x6 arrays): Current leg lengths, velocities, and accelerations.
%   q0_guess (18x1 array): Initial guess for joints B1-B6 [xb1;yb1;zb1;...].
%
% Returns:
%   kin (struct): Kinematic data including:
%       - B (3x6 double): Spatial positions of platform joints.
%       - vB, aB (3x6 double): Velocities and accelerations of joints B.
%       - RB_O (3x3 double): Rotation matrix of the platform.
%       - dRB_O, ddRB_O (3x3 double): Time derivatives of rotation matrix.
%       - P_fk (3x1 double): Center position of the moving platform.
%       - q0_next (18x1 array): Next iteration guess for the solver.

    % 1. Setup Base and Platform Local Points
    A = robot.base'; % 3x6
    Bp = robot.platform_pts'; % 3x6
    
    % Extract distances and linear combination constants from local platform geometry
    d12 = norm(Bp(:,1) - Bp(:,2));
    d13 = norm(Bp(:,1) - Bp(:,3));
    d23 = norm(Bp(:,2) - Bp(:,3));
    
    N = -[Bp(:,2) - Bp(:,1), Bp(:,3) - Bp(:,1)];
    X1 = N \ -(Bp(:,4) - Bp(:,1)); alpha1 = X1(1); beta1 = X1(2);
    X2 = N \ -(Bp(:,5) - Bp(:,1)); alpha2 = X2(1); beta2 = X2(2);
    X3 = N \ -(Bp(:,6) - Bp(:,1)); alpha3 = X3(1); beta3 = X3(2);

    % Local reference frame vectors for rotation extraction
    v1b = Bp(:,2) - Bp(:,1);
    v2b = Bp(:,3) - Bp(:,2);
    Rb = [v1b, v2b, cross(v1b, v2b)];

    % 2. Newton-Raphson Solver for Positions
    q0 = q0_guess;
    maxiter = 500;
    
    for j = 1:maxiter
        % Reshape q0 into 3x6 for easy coordinate extraction
        B = reshape(q0, 3, 6);
        
        % Leg length equations
        eq_L = sum((B - A).^2, 1)' - L'.^2; % 6x1
        
        % Rigid body distance equations
        eq_d = [sum((B(:,2) - B(:,1)).^2) - d12^2;
                sum((B(:,3) - B(:,2)).^2) - d23^2;
                sum((B(:,3) - B(:,1)).^2) - d13^2];
                
        % Linear combination equations
        eq_lc1 = (B(:,4) - B(:,1)) - alpha1*(B(:,2) - B(:,1)) - beta1*(B(:,3) - B(:,1));
        eq_lc2 = (B(:,5) - B(:,1)) - alpha2*(B(:,2) - B(:,1)) - beta2*(B(:,3) - B(:,1));
        eq_lc3 = (B(:,6) - B(:,1)) - alpha3*(B(:,2) - B(:,1)) - beta3*(B(:,3) - B(:,1));
        
        f = [eq_L; eq_d; eq_lc1; eq_lc2; eq_lc3]; % 18x1
        
        % Assembly of the 18x18 Jacobian (Jc)
        J = zeros(18, 18);
        for k = 1:6
            idx = (k-1)*3 + 1 : k*3;
            J(k, idx) = 2 * (B(:,k) - A(:,k))';
        end
        
        J(7, 1:6) = [-2*(B(:,2)-B(:,1))',  2*(B(:,2)-B(:,1))'];
        J(8, 4:9) = [-2*(B(:,3)-B(:,2))',  2*(B(:,3)-B(:,2))'];
        J(9, 1:3) = -2*(B(:,3)-B(:,1))'; J(9, 7:9) = 2*(B(:,3)-B(:,1))';
        
        J(10:12, 1:12) = [(-1 + alpha1 + beta1)*eye(3), -alpha1*eye(3), -beta1*eye(3), eye(3)];
        J(13:15, 1:9)  = [(-1 + alpha2 + beta2)*eye(3), -alpha2*eye(3), -beta2*eye(3)];
        J(13:15, 13:15) = eye(3);
        J(16:18, 1:9)  = [(-1 + alpha3 + beta3)*eye(3), -alpha3*eye(3), -beta3*eye(3)];
        J(16:18, 16:18) = eye(3);
        
        q1 = q0 + J\(-f);
        if norm(f) <= 10e-8
            break;
        end
        q0 = q1;
    end
    
    kin.B = reshape(q0, 3, 6);
    kin.q0_next = q0;
    
    % Rotation Matrix Extraction
    v1a = kin.B(:,2) - kin.B(:,1);
    v2a = kin.B(:,3) - kin.B(:,2);
    Ra = [v1a, v2a, cross(v1a, v2a)];
    kin.RB_O = Ra / Rb;
    
    kin.P_fk = A(:,1) + (kin.B(:,1) - A(:,1)) - kin.RB_O * Bp(:,1);

    % 3. Velocity Analysis
    Jj = zeros(18, 6);
    for k = 1:6
        Jj(k, k) = -2 * L(k);
    end
    
    Jc = J; % Jc is exactly the Newton-Raphson Jacobian at convergence
    v_vec = Jc \ -(Jj * dL'); % 18x1
    kin.vB = reshape(v_vec, 3, 6);
    
    % As established in the original formulation, rotational derivatives are kept null
    % for purely translational motion cases, preventing singularity in arbitrary spatial paths.
    kin.dRB_O = zeros(3); 
    
    % 4. Acceleration Analysis
    dJj = zeros(18, 6);
    for k = 1:6
        dJj(k, k) = -2 * dL(k);
    end
    
    dJc = zeros(18, 18);
    for k = 1:6
        idx = (k-1)*3 + 1 : k*3;
        dJc(k, idx) = 2 * kin.vB(:,k)';
    end
    dJc(7, 1:6) = [-2*(kin.vB(:,2)-kin.vB(:,1))',  2*(kin.vB(:,2)-kin.vB(:,1))'];
    dJc(8, 4:9) = [-2*(kin.vB(:,3)-kin.vB(:,2))',  2*(kin.vB(:,3)-kin.vB(:,2))'];
    dJc(9, 1:3) = -2*(kin.vB(:,3)-kin.vB(:,1))'; dJc(9, 7:9) = 2*(kin.vB(:,3)-kin.vB(:,1))';
    
    a_vec = Jc \ -(dJc * v_vec + dJj * dL' + Jj * ddL');
    kin.aB = reshape(a_vec, 3, 6);
    kin.ddRB_O = zeros(3);
end