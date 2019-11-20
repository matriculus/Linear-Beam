clear;clc;
syms x L q m D A pA e31 b t e33 z real
syms w1 t1 w2 t2 real
syms x1 x2 real
syms p_inf U_inf M_inf

transverse = [w1,t1,w2,t2]';

Pi = [1,x];
Pmat = [subs(Pi,x1);subs(Pi,x2)];
Pn = (Pi/Pmat)';

P = Pn;

Ni = [1,x,x.^2,x.^3];
dNi = diff(Ni,x);
Nmat = [subs(Ni,x1);subs(dNi,x1);subs(Ni,x2);subs(dNi,x2)];
Nn = (Ni/Nmat)';

N = Nn;

dN = diff(N,x);
dP = diff(P,x);
ddN = diff(dN,x);

dwx0 = dN'*transverse;

K11 = int(A*(dP*dP'),x,x1,x2);
K12 = int(0.5*A*(dN'*transverse)*dP*dN',x,x1,x2);
K21 = int(A*dwx0*dN*dP',x,x1,x2);
K22 = int(D*(ddN*ddN') + (A*dwx0.^2)*(dN*dN'),x,x1,x2);
M11 = int(pA*(P*P'),x,x1,x2);
M22 = int(pA*N*N',x,x1,x2);

K = int(D*(ddN*ddN'),x,x1,x2);
Kfl = int(-b*p_inf*U_inf.^2./sqrt(M_inf.^2-1).*dN*N',x,x1,x2);

F1 = int([P,dP],x,x1,x2);
F21 = int(N,x,x1,x2);
F22 = int(dN,x,x1,x2);
F23 = int(ddN,x,x1,x2);
F2 = [F21,F22,F23];

P1 = int(e31*b*dP,x,x1,x2);
P2 = int(-e31*b*z*ddN,x,x1,x2);

Q1 = int(-t*e31/e33*dP',x,x1,x2);
Q2 = int(t*e31/e33*z*ddN',x,x1,x2);


NN = matlabFunction(N);
PP = matlabFunction(P);
FF = matlabFunction(F2);


%%
% plot([0:0.01:1]',NN(1,0:0.01:1)','LineWidth',2);
% xlabel('X axis')
% ylabel('Displacement');
% legend('N_1','N_2','N_3','N_4','Location','east');
% saveas(gcf,'shapes_w','fig');
% saveas(gcf,'shapes_w','png');
% 
% %%
% plot([0:0.01:1]',PP(1,0:0.01:1)','LineWidth',2);
% xlabel('X axis')
% ylabel('Displacement');
% legend('P_1','P_2','Location','east');
% saveas(gcf,'shapes_u','fig');
% saveas(gcf,'shapes_u','png');