clear;clc;
syms x L q m D A pA e31 b z real
syms w1 t1 w2 t2 real

transverse = [w1,t1,w2,t2]';

P = [1-x./L;x./L];
N = 1./L^3.*[2.*x.^3 - 3.*x.^2.*L + L^3;
    x.^3.*L-2.*x.^2.*L.^2+x.*L.^3;
    -2.*x.^3+3.*x.^2.*L;
    x.^3.*L-x.^2.*L.^2];

dN = diff(N,x);
dP = diff(P,x);
ddN = diff(dN,x);

dwx0 = dN'*transverse;

K11 = int(A*(dP*dP'),x,0,L);
K12 = int(0.5*A*(dN'*transverse)*dP*dN',x,0,L);
K21 = int(A*dwx0*dN*dP',x,0,L);
K22 = int(D*(ddN*ddN') + (A*dwx0.^2)*(dN*dN'),x,0,L);
M11 = int(pA*(P*P'),x,0,L);
M22 = int(pA*N*N',x,0,L);

F1 = int([P,dP],x,0,L);
F2 = int([N,dN,ddN],x,0,L);

P1 = int(e31*b*dP,x,0,L);
P2 = int(-e31*b*z*ddN,x,0,L);


NN = matlabFunction(N);
PP = matlabFunction(P);

%%
plot([0:0.01:1]',NN(1,0:0.01:1)','LineWidth',2);
xlabel('X axis')
ylabel('Displacement');
legend('N_1','N_2','N_3','N_4','Location','east');
saveas(gcf,'shapes_w','fig');
saveas(gcf,'shapes_w','png');

%%
plot([0:0.01:1]',PP(1,0:0.01:1)','LineWidth',2);
xlabel('X axis')
ylabel('Displacement');
legend('P_1','P_2','Location','east');
saveas(gcf,'shapes_u','fig');
saveas(gcf,'shapes_u','png');