function obj = get_actuator_matrices(material, element_coordinates, u, w, tz, voltage)

order = [1,3,4,2,5,6];
transverse = transpose([w(1),tz(1),w(2),tz(2)]);
gauss = get_gaussian(element_coordinates(:), 2);

J = [0,1]*(gauss.matrix\element_coordinates(:));

rho = material.rho;
A = material.A;
B = material.B;
D = material.D;
e31 = material.piezoelectric_constant;
width = material.w;
ra = material.lever_arm;

K11 = 0;
K12 = 0;
K21 = 0;
K22 = 0;

T22 = 0;

F1 = 0;
F2 = 0;

shape_functions = get_shape_functions(gauss.coordinates, element_coordinates);

for g=1:gauss.gp
    weight = gauss.weights(g);
    
    phi = shape_functions.w(g,:);
    dphi_x = shape_functions.dwx(g,:);
    ddphi_xx = shape_functions.ddwxx(g,:);
    zhi = shape_functions.u(g,:);
    dzhi_x = shape_functions.dux(g,:);
    
    w0 = phi*transverse;
    dwx0 = dphi_x*transverse;
    u0 = zhi*u;
    dux0 = dzhi_x*u;
    
    gx = gauss.coordinates(g);
    
    Jacobian = J * weight;
    
    dK11 = A*(dzhi_x'*dzhi_x);
    K11 = K11 + dK11 * Jacobian;
    
    dK12 = 0.5*(A*dwx0)*(dzhi_x'*dphi_x);
    K12 = K12 + dK12 * Jacobian;
    
    dK21 = (A*dwx0)*(dphi_x'*dzhi_x);
    K21 = K21 + dK21 * Jacobian;
    
    dK22 = D*(ddphi_xx'*ddphi_xx) + (A*dwx0^2)*(dphi_x'*dphi_x);
    K22 = K22 + dK22 * Jacobian;
    
    F1 = F1 + e31*width*zhi'*Jacobian;
    F2 = F2 - e31*width*ra*ddphi_xx'*Jacobian;
    
    dT22 = A*(dux0 + dwx0*dwx0)*(dphi_x'*dphi_x);
    T22 = T22 + (dK22 + dT22) * Jacobian;
    
end
T11 = K11;
T21 = K21;
T12 = 2*K12;

T = [T11, T12; T21, T22];
obj.tangent = T(order, order);
external_force = [F1;F2]*voltage;
internal_force = [K11,K12;K21,K22]*[u;transverse];
obj.residue = external_force(order) - internal_force(order);
end