function torques = compute_inverse_dynamics_6sps(robot, kin, masses)
% COMPUTE_INVERSE_DYNAMICS_6SPS Calculates required actuation forces.
%
% Applies D'Alembert's principle with 3D equimomental tetrahedrons and 
% Davies' method (graph theory cutsets) to compute required actuation forces 
% for the 6 sliding legs (SPS joints) of the Stewart-Gough platform.
%
% Args:
%   robot (struct): Structural and geometric parameters.
%   kin (struct): Kinematic states of the 6 legs and platform.
%   masses (struct): Kinematics of the 13 equimomental tetrahedrons.
%
% Returns:
%   torques (struct): Actuation forces for the 6 legs (F1 to F6 in Newtons).

    %=====================================================================%
    %               1. INERTIAL & WEIGHT SCREWS (13 BODIES)               %
    %=====================================================================%
    S_in = zeros(6, 13);
    S_w  = zeros(6, 13);
    mag_in = zeros(13, 1);
    mag_w  = zeros(13, 1);

    for i = 1:6
        % --- Part 1 (e.g., Cylinder) ---
        idx1 = 2*i - 1;
        m1 = robot.legs(i).part1.m;
        
        S_temp1 = zeros(6,1);
        for j = 1:4
            p = masses.part1.p_O(:, j, i);
            a = masses.part1.ap_O(:, j, i);
            S_temp1 = S_temp1 + [cross(p, a); a];
        end
        Sin_full1 = (-m1/4) * S_temp1;
        so1 = cross(Sin_full1(4:6), Sin_full1(1:3)) / norm(Sin_full1(4:6))^2;
        mag_in(idx1) = norm(Sin_full1(4:6));
        us1 = Sin_full1(4:6) / norm(Sin_full1(4:6));
        S_in(:, idx1) = [cross(so1, us1); us1];
        
        cm1 = mean(masses.part1.p_O(:, :, i), 2);
        mag_w(idx1) = m1 * robot.g;
        S_w(:, idx1) = [-cm1(2); cm1(1); 0; 0; 0; -1];

        % --- Part 2 (e.g., Piston) ---
        idx2 = 2*i;
        m2 = robot.legs(i).part2.m;
        
        S_temp2 = zeros(6,1);
        for j = 1:4
            p = masses.part2.p_O(:, j, i);
            a = masses.part2.ap_O(:, j, i);
            S_temp2 = S_temp2 + [cross(p, a); a];
        end
        Sin_full2 = (-m2/4) * S_temp2;
        so2 = cross(Sin_full2(4:6), Sin_full2(1:3)) / norm(Sin_full2(4:6))^2;
        mag_in(idx2) = norm(Sin_full2(4:6));
        us2 = Sin_full2(4:6) / norm(Sin_full2(4:6));
        S_in(:, idx2) = [cross(so2, us2); us2];
        
        cm2 = mean(masses.part2.p_O(:, :, i), 2);
        mag_w(idx2) = m2 * robot.g;
        S_w(:, idx2) = [-cm2(2); cm2(1); 0; 0; 0; -1];
    end

    % --- Moving Platform (Body 13) ---
    m_plat = robot.platform.m;
    S_temp_p = zeros(6,1);
    for j = 1:4
        p = masses.plat.p_O(:, j);
        a = masses.plat.ap_O(:, j);
        S_temp_p = S_temp_p + [cross(p, a); a];
    end
    Sin_full_p = (-m_plat/4) * S_temp_p;
    so_p = cross(Sin_full_p(4:6), Sin_full_p(1:3)) / norm(Sin_full_p(4:6))^2;
    mag_in(13) = norm(Sin_full_p(4:6));
    us_p = Sin_full_p(4:6) / norm(Sin_full_p(4:6));
    S_in(:, 13) = [cross(so_p, us_p); us_p];
    
    cm_plat = mean(masses.plat.p_O, 2);
    mag_w(13) = m_plat * robot.g;
    S_w(:, 13) = [-cm_plat(2); cm_plat(1); 0; 0; 0; -1];

    %=====================================================================%
    %               2. PRIMARY VARIABLES & JOINT SCREWS                   %
    %=====================================================================%
    s_A = zeros(6, 18);
    s_B = zeros(6, 18);
    s_PT = zeros(6, 36); % Interleaved Prismatic Forces and Torques
    
    for i = 1:6
        % Base Joints (A1-A6)
        A = robot.base(i, :)';
        idxA = (i-1)*3 + 1;
        s_A(:, idxA)   = [0; A(3); -A(2); 1; 0; 0];
        s_A(:, idxA+1) = [-A(3); 0; A(1); 0; 1; 0];
        s_A(:, idxA+2) = [A(2); -A(1); 0; 0; 0; 1];

        % Platform Joints (B1-B6)
        B = kin.B(:, i);
        s_B(:, idxA)   = [0; B(3); -B(2); 1; 0; 0];
        s_B(:, idxA+1) = [-B(3); 0; B(1); 0; 1; 0];
        s_B(:, idxA+2) = [B(2); -B(1); 0; 0; 0; 1];

        % Prismatic Actuation and Restriction Screws (P1-P6)
        Lv = B - A;
        s_vec = Lv / norm(Lv);
        cG = s_vec(3); sG = sqrt(s_vec(1)^2 + s_vec(2)^2);
        if sG == 0, sBt = 0; cBt = 1; else, sBt = s_vec(2)/sG; cBt = s_vec(1)/sG; end
        
        RA_O = [cBt*cG, -sBt, cBt*sG; 
                sBt*cG,  cBt, sBt*sG; 
               -sG,      0,   cG];
               
        So = A + s_vec;
        Fpx = RA_O*[1;0;0]; Fpy = RA_O*[0;1;0]; Fpz = RA_O*[0;0;1];
        
        idxP = (i-1)*6 + 1;
        s_PT(:, idxP)   = [cross(So, Fpx); Fpx];
        s_PT(:, idxP+1) = [cross(So, Fpy); Fpy];
        s_PT(:, idxP+2) = [cross(So, Fpz); Fpz]; % Actuation Force Direction
        s_PT(:, idxP+3) = [Fpx; 0;0;0];
        s_PT(:, idxP+4) = [Fpy; 0;0;0];
        s_PT(:, idxP+5) = [Fpz; 0;0;0];
    end

    % Flat assembly matching the exact column order of the Q Cutset matrix
    A_d = [s_A, s_B, s_PT, S_in, S_w];

    %=====================================================================%
    %               3. DAVIES METHOD (GRAPH CUTSET MATRIX)                %
    %=====================================================================%
    
    Q1=[   1      1      1      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      1      1      1      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      1      1      1      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      1      1      1      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      1      1      1      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      1      1      1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0               ]    ; 
   
    Q2=[   0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           1      1      1      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      1      1      1      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      1      1      1      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      1      1      1      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      1      1      1      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      1      1      1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0               ]    ; 
   
    Q3=[   0      0      0      0      0      0      0      0      0      0      0      0      1      1      1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0     -1     -1     -1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      1      1      1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0     -1     -1     -1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           1      1      1      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      1      1      1      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      1      1      1      1      1      1      1      1      1               ]    ; 
       
    Q4=[   1      1      1      1      1      1      1      1      1      1      1      1      1      1      1                    ; 
          -1     -1     -1      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0     -1     -1     -1     -1     -1     -1      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0     -1     -1     -1     -1     -1     -1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           1      1      1      1      1      1      1      1      1      1      1      1      1      1      1                    ; 
          -1     -1     -1      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0     -1     -1     -1     -1     -1     -1      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0     -1     -1     -1     -1     -1     -1                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           0      0      0      0      0      0      0      0      0      0      0      0      0      0      0                    ; 
           1      1      1      1      1      1      1      1      1      1      1      1      1      1      1               ]    ; 
   
    Q5=[   1      1      1      1      1      1      1      1      1      1      1      1       1      1      0                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      1                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      0                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      0                  ; 
          -1     -1     -1     -1     -1     -1      0      0      0      0      0      0       0      0      0                  ; 
           0      0      0      0      0      0     -1     -1     -1     -1     -1     -1       0      0      0                  ; 
           1      1      1      1      1      1      1      1      1      1      1      1       0      0      0                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      0                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      0                  ; 
           0      0      0      0      0      0      0      0      0      0      0      0       0      0      0                  ; 
          -1     -1     -1     -1     -1     -1      0      0      0      0      0      0       0      0      0                  ; 
           0      0      0      0      0      0     -1     -1     -1     -1     -1     -1       0      0      0                  ; 
           1      1      1      1      1      1      1      1      1      1      1      1       0      1      0              ]    ; 
   
    Q6=[    1       0       1       0       1       0       1       0       1       1      1      1       0      1      0        ; 
            0       0       0       0       0       0       0       0       0       0      0      0       1      0      0        ; 
            0       1       0       0       0       0       0       0       0       0      0      0       0      0      1        ; 
            0       0       0       1       0       0       0       0       0       0      0      0       0      0      0        ; 
            0       0       0       0       0       1       0       0       0       0      0      0       0      0      0        ; 
            0       0       0       0       0       0       0       1       0       0      0      0       0      0      0        ; 
            1       0       1       0       1       0       1       0       1       1      0      0       0      1      0        ; 
           -1       0       0       0       0       0       0       0       0       0      0      0       0     -1      0        ; 
            0       0      -1       0       0       0       0       0       0       0      0      0       0      0      0        ; 
            0       0       0       0      -1       0       0       0       0       0      0      0       0      0      0        ; 
            0       0       0       0       0       0      -1       0       0       0      0      0       0      0      0        ; 
            0       0       0       0       0       0       0       0      -1       0      0      0       0      0      0        ; 
            1       0       1       0       1       0       1       0       1       1      0      1       0      1      0   ]    ; 
    
    Q7=[   1      0       1      0      1      0      1      1                                                                 ; 
           0      0       0      0      0      0      0      0                                                                 ; 
           0      0       0      0      0      0      0      0                                                                 ; 
           0      1       0      0      0      0      0      0                                                                 ; 
           0      0       0      1      0      0      0      0                                                                 ; 
           0      0       0      0      0      1      0      0                                                                 ; 
           1      0       1      0      1      0      1      1                                                                 ; 
           0      0       0      0      0      0      0      0                                                                 ; 
          -1      0       0      0      0      0      0      0                                                                 ; 
           0      0      -1      0      0      0      0      0                                                                 ; 
           0      0       0      0     -1      0      0      0                                                                 ; 
           0      0       0      0      0      0     -1      0                                                                 ; 
           1      0       1      0      1      0      1      1                                                              ]    ;   

    Q = [Q1 Q2 Q3 Q4 Q5 Q6 Q7];

    % Network Actions Matrix Allocation
    A_n = [A_d*diag(Q(1,:));
           A_d*diag(Q(2,:));
           A_d*diag(Q(3,:));
           A_d*diag(Q(4,:));
           A_d*diag(Q(5,:));
           A_d*diag(Q(6,:));
           A_d*diag(Q(7,:));
           A_d*diag(Q(8,:));
           A_d*diag(Q(9,:));
           A_d*diag(Q(10,:));
           A_d*diag(Q(11,:));
           A_d*diag(Q(12,:));
           A_d*diag(Q(13,:))];

    A_n_s = A_n(:, 1:72);   % Secondary Variables (Unknown Constraints)
    A_n_p = A_n(:, 73:end); % Primary Variables (Known D'Alembert & Weights)
    
    % Magnitudes array flattened
    psi_p = [mag_in; mag_w]; 

    % Linear solver for the unknown constraints and torques
    Psi = A_n_s \ (-A_n_p * psi_p);
    
    % The required actuation forces correspond to the local z-axis (Fpz) 
    % of the prismatic joints in the secondary variables array.
    % Indices: Leg 1 (39), Leg 2 (45), Leg 3 (51), Leg 4 (57), Leg 5 (63), Leg 6 (69)
    torques.F1 = Psi(39);
    torques.F2 = Psi(45);
    torques.F3 = Psi(51);
    torques.F4 = Psi(57);
    torques.F5 = Psi(63);
    torques.F6 = Psi(69);
end