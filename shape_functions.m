syms x len q m real

N = 1./len^3.*[2.*x.^3 - 3.*x.^2.*len + len^3;
    x.^3.*len-2.*x.^2.*len.^2+x.*len.^3;
    -2.*x.^3+3.*x.^2.*len;
    x.^3.*len-x.^2.*len.^2];

dN = diff(N,x);
ddN = diff(dN,x);

F = int([N,dN,ddN],x,0,len);