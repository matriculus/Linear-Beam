function obj = get_shape_functions(polynomial, matrices)
transverse = matrices.transverse;
axial = matrices.axial;

obj.w = polynomial.w/transverse;
obj.dwx = polynomial.dwx/transverse;
obj.ddwxx = polynomial.ddwxx/transverse;
obj.u = polynomial.u/axial;
obj.dux = polynomial.dux/axial;
end