function obj = get_actuator_matrices(EI, rho, e31b, len, ra)

obj.stiffness = EI/len^3*[12, 6*len, -12, 6*len;
    6*len, 4*len^2, -6*len, 2*len^2;
    -12, -6*len, 12, -6*len;
    6*len, 2*len^2, -6*len, 4*len^2];

obj.mass = rho*len/420*[156, 22*len, 54, -13*len;
    22*len, 4*len^2, 13*len, -3*len^2;
    54, 13*len, 156, -22*len;
    -13*len, -3*len^2, -22*len, 4*len^2];

obj.force = [0; e31b*ra; 0; -e31b*ra];
end