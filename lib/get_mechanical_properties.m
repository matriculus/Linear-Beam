function obj = get_mechanical_properties(E, nu, rho, w, d1, d2)
E1 = @(x)x.*0 + 1;
E2 = @(x)x;
E3 = @(x)x.*x;
I = integral(E3,d1,d2);
obj.dA = abs(w)*abs((d1-d2));
obj.A = E * integral(E1,d1,d2)*w;
obj.B = E * integral(E2,d1,d2)*w;
if abs(obj.B) < 1e-8
    obj.B = 0;
end
obj.D = E*I*w;
obj.w = w;
obj.rho = rho*obj.dA;
obj.lever_arm = 0.5*(d2-d1)+d1;
end