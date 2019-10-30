function obj = get_shape_matrices(x)
ply = get_polynomial(x);
obj.axial = ply.u;
obj.transverse = [ply.w(1,:);ply.dwx(1,:);ply.w(2,:);ply.dwx(2,:)];
end