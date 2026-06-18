function kin = compute_kinematics_3rrr(robot, Px, Py, dth, ddth, q0_guess)
% COMPUTE_KINEMATICS_3RRR Forward/Inverse kinematics for the 3-RRR robot.
%
% Uses analytical solutions for the active joints and a Newton-Raphson
% numerical solver for the passive joint closures.
%
% Args:
%   robot (struct): Structural parameters.
%   Px, Py (double): End-effector target position.
%   dth, ddth (3x1 array): Joint velocities and accelerations for active joints.
%   q0_guess (6x1 array): Initial guess for the Newton-Raphson solver [xc1;yc1;xc2;yc2;xc3;yc3].
%
% Returns:
%   kin (struct): Positions, velocities, and accelerations of all joints, 
%                 along with the updated guess for the next iteration.

    % Extract base points for readability
    xa1 = robot.base(1,1); ya1 = robot.base(1,2);
    xa2 = robot.base(2,1); ya2 = robot.base(2,2);
    xa3 = robot.base(3,1); ya3 = robot.base(3,2);
    
    d1 = robot.d(1); d2 = robot.d(2); d3 = robot.d(3);
    he = (((3)^(0.5))/2)*d1;
    phi = 0; % Platform rotation

    % 1. Inverse Kinematics (Find Active Joint Angles th1, th2, th3)
    xC1 = Px - d1/2; 
    yC1 = Py - he/3;
    
    % Leg 1 Inverse
    e11 = -2*yC1*robot.L(1); e12 = -2*xC1*robot.L(1); e13 = xC1^2 + yC1^2 + robot.L(1)^2 - robot.L(2)^2;
    r1 = roots([e13-e12, 2*e11, e13+e12]);
    th1 = 2*atan(r1(1));
    
    % Leg 2 Inverse
    e21 = -2*yC1*robot.L(3) + 2*ya2*robot.L(3) - 2*robot.L(3)*d1*sin(phi);
    e22 = -2*xC1*robot.L(3) - 2*robot.L(3)*d1*cos(phi) + 2*xa2*robot.L(3);
    e23 = xC1^2 + yC1^2 - 2*xC1*xa2 - 2*yC1*ya2 + xa2^2 + ya2^2 + d1^2 + robot.L(3)^2 - robot.L(4)^2 ...
          + 2*xC1*d1*cos(phi) + 2*yC1*d1*sin(phi) - 2*xa2*d1*cos(phi) - 2*ya2*d1*sin(phi);
    r2 = roots([e23-e22, 2*e21, e23+e22]);
    th2 = 2*atan(r2(2));
    
    % Leg 3 Inverse
    e31 = -2*yC1*robot.L(5) + 2*ya3*robot.L(5) - 2*robot.L(5)*d1*sin(phi+(pi/3));
    e32 = -2*xC1*robot.L(5) - 2*robot.L(5)*d1*cos(phi+(pi/3)) + 2*xa3*robot.L(5);
    e33 = xC1^2 + yC1^2 - 2*xC1*xa3 - 2*yC1*ya3 + xa3^2 + ya3^2 + d1^2 + robot.L(5)^2 - robot.L(6)^2 ...
          + 2*xC1*d1*cos(phi+(pi/3)) + 2*yC1*d1*sin(phi+(pi/3)) - 2*xa3*d1*cos(phi+(pi/3)) - 2*ya3*d1*sin(phi+(pi/3));
    r3 = roots([e33-e32, 2*e31, e33+e32]);
    th3 = 2*atan(r3(2));
    
    % Active joint position vectors (B1, B2, B3)
    xb1 = xa1 + robot.L(1)*cos(th1); yb1 = ya1 + robot.L(1)*sin(th1);
    xb2 = xa2 + robot.L(3)*cos(th2); yb2 = ya2 + robot.L(3)*sin(th2);
    xb3 = xa3 + robot.L(5)*cos(th3); yb3 = ya3 + robot.L(5)*sin(th3);

    % 2. Forward Kinematics (Newton-Raphson Solver for C1, C2, C3)
    q0 = q0_guess;
    maxiter = 500;
    
    for j = 1:maxiter
        xc1 = q0(1); yc1 = q0(2);
        xc2 = q0(3); yc2 = q0(4);
        xc3 = q0(5); yc3 = q0(6);
        
        eq1 = (xc1-xb1)^2 + (yc1-yb1)^2 - robot.L(2)^2;
        eq2 = (xc2-xb2)^2 + (yc2-yb2)^2 - robot.L(4)^2;
        eq3 = (xc3-xb3)^2 + (yc3-yb3)^2 - robot.L(6)^2;
        eq4 = (xc2-xc1)^2 + (yc2-yc1)^2 - d1^2;
        eq5 = (xc3-xc2)^2 + (yc3-yc2)^2 - d2^2;
        eq6 = (xc1-xc3)^2 + (yc1-yc3)^2 - d3^2;
        
        f = [eq1; eq2; eq3; eq4; eq5; eq6];
        
        J_NR = [ 2*(xc1-xb1), 2*(yc1-yb1), 0, 0, 0, 0;
                 0, 0, 2*(xc2-xb2), 2*(yc2-yb2), 0, 0;
                 0, 0, 0, 0, 2*(xc3-xb3), 2*(yc3-yb3);
                -2*(xc2-xc1), -2*(yc2-yc1), 2*(xc2-xc1), 2*(yc2-yc1), 0, 0;
                 0, 0, -2*(xc3-xc2), -2*(yc3-yc2), 2*(xc3-xc2), 2*(yc3-yc2);
                 2*(xc1-xc3), 2*(yc1-yc3), 0, 0, -2*(xc1-xc3), -2*(yc1-yc3)];
                 
        q1 = q0 + J_NR\(-f);
        if norm(f) <= 10e-8
            break;
        end
        q0 = q1;
    end
    
    % Store Positions
    kin.th = [th1; th2; th3];
    kin.B = [xb1, xb2, xb3; yb1, yb2, yb3];
    kin.C = [q0(1), q0(3), q0(5); q0(2), q0(4), q0(6)];
    kin.P = [q0(1) + d1/2; q0(2) + he/3];
    kin.q0_next = q0; % Return for next iteration

    % 3. Differential Kinematics (Jacobians)
    Jj = [robot.L(1)*sin(th1), 0, 0;
         -robot.L(1)*cos(th1), 0, 0;
          0, robot.L(3)*sin(th2), 0;
          0,-robot.L(3)*cos(th2), 0;
          0, 0, robot.L(5)*sin(th3);
          0, 0,-robot.L(5)*cos(th3);
          zeros(6,3)];
          
    Ja = [-2*(q0(1)-xb1), -2*(q0(2)-yb1), 0, 0, 0, 0, 2*(q0(1)-xb1), 2*(q0(2)-yb1), 0, 0, 0, 0;
           0, 0, -2*(q0(3)-xb2), -2*(q0(4)-yb2), 0, 0, 0, 0, 2*(q0(3)-xb2), 2*(q0(4)-yb2), 0, 0;
           0, 0, 0, 0, -2*(q0(5)-xb3), -2*(q0(6)-yb3), 0, 0, 0, 0, 2*(q0(5)-xb3), 2*(q0(6)-yb3);
           0, 0, 0, 0, 0, 0, -2*(q0(3)-q0(1)), -2*(q0(4)-q0(2)), 2*(q0(3)-q0(1)), 2*(q0(4)-q0(2)), 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, -2*(q0(5)-q0(3)), -2*(q0(6)-q0(4)), 2*(q0(5)-q0(3)), 2*(q0(6)-q0(4));
           0, 0, 0, 0, 0, 0, 2*(q0(1)-q0(5)), 2*(q0(2)-q0(6)), 0, 0, -2*(q0(1)-q0(5)), -2*(q0(2)-q0(6))];

    Jc = [eye(6), zeros(6,6); Ja];
    
    % Velocities
    v = Jc \ -(Jj * dth);
    kin.vB = [v(1), v(3), v(5); v(2), v(4), v(6)];
    kin.vC = [v(7), v(9), v(11); v(8), v(10), v(12)];
    kin.vP = [mean(kin.vC(1,:)); mean(kin.vC(2,:))];

    % Accelerations
    dJi = [robot.L(1)*cos(th1)*dth(1), 0, 0;
           robot.L(1)*sin(th1)*dth(1), 0, 0;
           0, robot.L(3)*cos(th2)*dth(2), 0;
           0, robot.L(3)*sin(th2)*dth(2), 0;
           0, 0, robot.L(5)*cos(th3)*dth(3);
           0, 0, robot.L(5)*sin(th3)*dth(3)];
           
    dJd = [-2*(kin.vC(1,1)-kin.vB(1,1)), -2*(kin.vC(2,1)-kin.vB(2,1)), 0, 0, 0, 0, 2*(kin.vC(1,1)-kin.vB(1,1)), 2*(kin.vC(2,1)-kin.vB(2,1)), 0, 0, 0, 0;
            0, 0, -2*(kin.vC(1,2)-kin.vB(1,2)), -2*(kin.vC(2,2)-kin.vB(2,2)), 0, 0, 0, 0, 2*(kin.vC(1,2)-kin.vB(1,2)), 2*(kin.vC(2,2)-kin.vB(2,2)), 0, 0;
            0, 0, 0, 0, -2*(kin.vC(1,3)-kin.vB(1,3)), -2*(kin.vC(2,3)-kin.vB(2,3)), 0, 0, 0, 0, 2*(kin.vC(1,3)-kin.vB(1,3)), 2*(kin.vC(2,3)-kin.vB(2,3));
            0, 0, 0, 0, 0, 0, -2*(kin.vC(1,2)-kin.vC(1,1)), -2*(kin.vC(2,2)-kin.vC(2,1)), 2*(kin.vC(1,2)-kin.vC(1,1)), 2*(kin.vC(2,2)-kin.vC(2,1)), 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, -2*(kin.vC(1,3)-kin.vC(1,2)), -2*(kin.vC(2,3)-kin.vC(2,2)), 2*(kin.vC(1,3)-kin.vC(1,2)), 2*(kin.vC(2,3)-kin.vC(2,2));
            0, 0, 0, 0, 0, 0, 2*(kin.vC(1,1)-kin.vC(1,3)), 2*(kin.vC(2,1)-kin.vC(2,3)), 0, 0, -2*(kin.vC(1,1)-kin.vC(1,3)), -2*(kin.vC(2,1)-kin.vC(2,3))];
            
    dJ = [zeros(6,12), dJi; dJd, zeros(6,3)];
    ac = Jc \ (-dJ * [v; dth] - Jj * ddth);
    
    kin.aB = [ac(1), ac(3), ac(5); ac(2), ac(4), ac(6)];
    kin.aC = [ac(7), ac(9), ac(11); ac(8), ac(10), ac(12)];
    kin.aP = [mean(kin.aC(1,:)); mean(kin.aC(2,:))];
end