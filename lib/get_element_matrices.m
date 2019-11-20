function obj = get_element_matrices(EI, rho, x1,x2)

Kfl = @(M_inf,U_inf,b,p_inf,x1,x2)[(U_inf^2*b*p_inf)/(2*(M_inf^2 - 1)^(1/2)), -(U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)),  (U_inf^2*b*p_inf)/(2*(M_inf^2 - 1)^(1/2)),   (U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2));
    (U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)), 0, -(U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)), -(b*p_inf*U_inf^2*x1^2 - 2*b*p_inf*U_inf^2*x1*x2 + b*p_inf*U_inf^2*x2^2)/(60*(M_inf^2 - 1)^(1/2));
    -(U_inf^2*b*p_inf)/(2*(M_inf^2 - 1)^(1/2)),  (U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)), -(U_inf^2*b*p_inf)/(2*(M_inf^2 - 1)^(1/2)),  -(U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2));
    -(U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)), (b*p_inf*U_inf^2*x1^2 - 2*b*p_inf*U_inf^2*x1*x2 + b*p_inf*U_inf^2*x2^2)/(60*(M_inf^2 - 1)^(1/2)),  (U_inf^2*b*p_inf*x1 - U_inf^2*b*p_inf*x2)/(10*(M_inf^2 - 1)^(1/2)),  0];

K = @(D,x1,x2)[ -(12*D)/(x1 - x2)^3,  (6*D)/(x1 - x2)^2,  (12*D)/(x1 - x2)^3,  (6*D)/(x1 - x2)^2;
    (6*D)/(x1 - x2)^2,   -(4*D)/(x1 - x2),  -(6*D)/(x1 - x2)^2,   -(2*D)/(x1 - x2);
    (12*D)/(x1 - x2)^3, -(6*D)/(x1 - x2)^2, -(12*D)/(x1 - x2)^3, -(6*D)/(x1 - x2)^2;
    (6*D)/(x1 - x2)^2,   -(2*D)/(x1 - x2),  -(6*D)/(x1 - x2)^2,   -(4*D)/(x1 - x2)];
M = @(pA,x1,x2)[-(13*pA*(x1 - x2))/35, (11*pA*(x1 - x2)^2)/210, -(9*pA*(x1 - x2))/70, -(13*pA*(x1 - x2)^2)/420;
    (11*pA*(x1 - x2)^2)/210,   -(pA*(x1 - x2)^3)/105,  (13*pA*(x1 - x2)^2)/420, (pA*(x1 - x2)^3)/140;
    -(9*pA*(x1 - x2))/70, (13*pA*(x1 - x2)^2)/420,    -(13*pA*(x1 - x2))/35, -(11*pA*(x1 - x2)^2)/210;
    -(13*pA*(x1 - x2)^2)/420,    (pA*(x1 - x2)^3)/140, -(11*pA*(x1 - x2)^2)/210,    -(pA*(x1 - x2)^3)/105];
F = @(x1,x2)[ x2/2 - x1/2, -1,  0;
    (x1 - x2)^2/12,  0, -1;
    x2/2 - x1/2,  1,  0;
    -(x1 - x2)^2/12,  0,  1];

obj.stiffness = K(EI,x1,x2) + 0.*Kfl(1.5,100,0.03,1,x1,x2);
obj.mass = M(rho,x1,x2);
obj.force = F(x1,x2);
end