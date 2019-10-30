function obj = get_mechanical_properties(E, nu, rho, w1, w2, d1, d2)
E1 = @(x)x.*0 + 1;
E2 = @(x)x;
E3 = @(x)x.*x;
b = integral(E1,w1,w2);
I = integral(E3,d1,d2);
obj.dA = abs((w1-w2))*abs((d1-d2));
obj.A = E * integral(E1,d1,d2)*b;
obj.B = E * integral(E2,d1,d2)*b;
obj.D = E*I*b;
obj.rho = rho*obj.dA;
end