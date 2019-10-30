function obj = get_polynomial(x)
O = zeros(size(x));
I = ones(size(x));
obj.w = [I,x,x.^2,x.^3];
obj.dwx = [O,I,2.*x,3.*x.^2];
obj.ddwxx = [O,O,2.*I,6.*x];
obj.u = [I,x];
obj.dux = [O,I];
end