function obj = get_shape_functions(values, coordinates)
matrices = get_shape_matrices(coordinates(:));
polynomial = get_polynomial(values);

transverse = matrices.transverse;
axial = matrices.axial;

obj.w = polynomial.w/transverse;
obj.dwx = polynomial.dwx/transverse;
obj.ddwxx = polynomial.ddwxx/transverse;
obj.u = polynomial.u/axial;
obj.dux = polynomial.dux/axial;
end