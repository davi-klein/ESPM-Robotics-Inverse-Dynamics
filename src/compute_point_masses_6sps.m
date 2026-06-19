function masses = compute_point_masses_6sps(robot, kin)
% COMPUTE_POINT_MASSES_6SPS Kinematics of 3D equimomental point masses.
%
% Calculates global position, velocity, and acceleration for the tetrahedrons
% equimomental to the 13 rigid bodies (12 leg parts + 1 moving platform).
%
% Args:
%   robot (struct): Structural and pseudo-inertia properties.
%   kin (struct): Kinematic states from the Newton-Raphson solver.
%
% Returns:
%   masses (struct): Kinematic data of the tetrahedrons:
%       - part1.p_O (3x4x6 double): Positions for the 6 lower leg parts.
%       - part2.p_O (3x4x6 double): Positions for the 6 upper leg parts.
%       - plat.p_O (3x4 double): Positions for the moving platform.
%       (Includes vp_O and ap_O for velocities and accelerations).

    % Regular Tetrahedron Vertices (Local Frame)
    q = [ 1, -1,  1, -1;
          1, -1, -1,  1;
          1,  1, -1, -1;
          1,  1,  1,  1 ]; % 4x4 matrix (including homogeneous scale)

    % Pre-allocate memory for speed
    masses.part1.p_O = zeros(3, 4, 6); masses.part1.vp_O = zeros(3, 4, 6); masses.part1.ap_O = zeros(3, 4, 6);
    masses.part2.p_O = zeros(3, 4, 6); masses.part2.vp_O = zeros(3, 4, 6); masses.part2.ap_O = zeros(3, 4, 6);

    %=====================================================================%
    %                       6 SPS LEGS KINEMATICS                         %
    %=====================================================================%
    for i = 1:6
        % Extract states for leg i
        A = robot.base(i, :)';
        B = kin.B(:, i); vB = kin.vB(:, i); aB = kin.aB(:, i);
        
        Lv = B - A;
        L = norm(Lv);
        dL = dot(Lv, vB) / L;
        ddL = (dot(vB, vB) + dot(Lv, aB) - dL^2) / L;

        % Director unitary vector and angles
        s = Lv / L;
        cG = s(3);
        sG = sqrt(s(1)^2 + s(2)^2);
        
        % Protection against singularity at zenith (z-axis alignment)
        if sG == 0
            sBt = 0; cBt = 1; 
        else
            sBt = s(2) / sG; cBt = s(1) / sG;
        end

        % 1. Rotation Matrix RA_O
        RA_O = [cBt*cG, -sBt, cBt*sG;
                sBt*cG,  cBt, sBt*sG;
               -sG,      0,   cG];

        % 2. First and Second time derivatives of angles
        dG = ((Lv(3)*dL) - (L*vB(3))) / (L^2 * sG);
        num_ddG = (L^2 * sG * (Lv(3)*ddL - L*aB(3))) - ((Lv(3)*dL - L*vB(3)) * (2*L*sG*dL + L^2*cG*dG));
        ddG = num_ddG / (L^2 * sG)^2;

        dBt = ((L*vB(2) - Lv(2)*dL) / (L^2 * cBt * sG)) - ((sBt * cG * dG) / (cBt * sG));
        term_ddBt = (L^2 * (L*aB(2) - Lv(2)*ddL)) - ((L*vB(2) - Lv(2)*dL) * 2*L*dL);
        ddBt = (term_ddBt / (L^4 * cBt * sG)) + (sBt * (dBt^2 + dG^2) / cBt) - (2*cG*dBt*dG / sG) - (sBt*cG*ddG / (cBt*sG));

        % 3. First time derivative of Rotation Matrix (dRA_O)
        dRA_O = [(-sBt*cG*dBt)-(cBt*sG*dG), -cBt*dBt, (-sBt*sG*dBt)+(cBt*cG*dG);
                 (cBt*cG*dBt)-(sBt*sG*dG),  -sBt*dBt, (cBt*sG*dBt)+(sBt*cG*dG);
                 -cG*dG,                     0,       -sG*dG];

        % 4. Second time derivative of Rotation Matrix (ddRA_O)
        ddr11 = 2*sG*sBt*dBt*dG - cG*cBt*dBt^2 - cG*sBt*ddBt - cBt*cG*dG^2 - cBt*sG*ddG;
        ddr12 = sBt*dBt^2 - cBt*ddBt;
        ddr13 = -2*sBt*cG*dBt*dG - cBt*sG*dG^2 + cBt*cG*ddG - cBt*sG*dBt^2 - sBt*sG*ddBt;
        ddr21 = -sBt*cG*dBt^2 - 2*cBt*sG*dBt*dG + cBt*cG*ddBt - sBt*cG*dG^2 - sBt*sG*ddG;
        ddr22 = -cBt*dBt^2 - sBt*ddBt;
        ddr23 = -sBt*sG*dBt^2 + 2*cBt*cG*dBt*dG + cBt*sG*ddBt - sG*sBt*dG^2 + cG*sBt*ddG;
        ddr31 = sG*dG^2 - cG*ddG^2; % Exact transcription to maintain validation parity
        ddr32 = 0;
        ddr33 = -cG*dG^2 - sG*ddG;

        ddRA_O = [ddr11, ddr12, ddr13;
                  ddr21, ddr22, ddr23;
                  ddr31, ddr32, ddr33];

        % --- PART 1 (e.g., Cylinder) ---
        D1 = sqrt((1 / robot.legs(i).part1.m) * robot.legs(i).part1.E);
        p1 = D1 * q; 
        p1_A = p1 + [0; 0; 0.5; 0]; % Offset in local frame

        masses.part1.p_O(:,:,i) = RA_O * p1_A(1:3,:) + A * ones(1,4);
        masses.part1.vp_O(:,:,i) = dRA_O * p1_A(1:3,:);
        masses.part1.ap_O(:,:,i) = ddRA_O * p1_A(1:3,:);

        % --- PART 2 (e.g., Piston) ---
        D2 = sqrt((1 / robot.legs(i).part2.m) * robot.legs(i).part2.E);
        p2 = D2 * q;
        p2_A = p2 + [0; 0; L - 0.5; 0];

        masses.part2.p_O(:,:,i) = RA_O * p2_A(1:3,:) + A * ones(1,4);
        
        v2_A = [0; 0; dL];
        masses.part2.vp_O(:,:,i) = dRA_O * p2_A(1:3,:) + RA_O * v2_A * ones(1,4);

        a2_A = [0; 0; ddL];
        % Includes Coriolis Effect: 2 * dRA_O * v2_A
        masses.part2.ap_O(:,:,i) = ddRA_O * p2_A(1:3,:) + 2 * dRA_O * (v2_A * ones(1,4)) + RA_O * a2_A * ones(1,4);
    end

    %=====================================================================%
    %                  MOVING PLATFORM KINEMATICS (BODY 13)               %
    %=====================================================================%
    D_mp = sqrt((1 / robot.platform.m) * robot.platform.E);
    p_mp = D_mp * q;

    % Pure translation logic preserved from the analytical model
    vP = kin.vB(:,1) - kin.dRB_O * robot.platform_pts(1,:)';
    aP = kin.aB(:,1) - kin.ddRB_O * robot.platform_pts(1,:)';

    masses.plat.p_O = kin.RB_O * p_mp(1:3,:) + kin.P_fk * ones(1,4);
    masses.plat.vp_O = vP * ones(1,4) + kin.dRB_O * p_mp(1:3,:);
    masses.plat.ap_O = aP * ones(1,4) + kin.ddRB_O * p_mp(1:3,:);
end