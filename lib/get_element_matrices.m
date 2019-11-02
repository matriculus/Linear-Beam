function obj = get_element_matrices(EI, rho, len)

obj.stiffness = EI/len^3*[12, 6*len, -12, 6*len;
    6*len, 4*len^2, -6*len, 2*len^2;
    -12, -6*len, 12, -6*len;
    6*len, 2*len^2, -6*len, 4*len^2];
obj.mass = rho*len/420*[156, 22*len, 54, -13*len;
    22*len, 4*len^2, 13*len, -3*len^2;
    54, 13*len, 156, -22*len;
    -13*len, -3*len^2, -22*len, 4*len^2];
obj.force = [0.5*len, -1, 0;1/12*len^2, 0, -1;0.5*len, 1, 0;-1/12*len^2, 0, 1];
end