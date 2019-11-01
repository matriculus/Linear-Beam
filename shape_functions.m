syms x len e31 ra w real

N = 1./len^3.*[2.*x.^3 - 3.*x.^2.*len + len^3;
    x.^3.*len-2.*x.^2.*len.^2+x.*len.^3;
    -2.*x.^3+3.*x.^2.*len;
    x.^3.*len-x.^2.*len.^2];

ddN = diff(diff(N,x),x);

F = int(-ra*e31*w*ddN,x,0,len);