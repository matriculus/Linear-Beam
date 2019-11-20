close all;
clear;clc;

n = 1;

kI = 150;
k1 = 100;
k2 = 10;
for i=1:n
    A = [0,1,0;0,0,1;-kI,-k1,-k2];
    B = [0,0,1]';
    C = [0,1,0];
    Q = eye(3);
    R = 1;
    [K,~,p] = lqr(A,B,Q,R);
    Ac = A-B*K;
    
    KK = -Ac(3,:);
    kI = KK(1);
    k1 = KK(2);
    k2 = KK(3);
end
disp(Ac);
disp(p);
fprintf('kI = %f;\nk1 = %f;\nk2 = %f;\n',KK(1),KK(2),KK(3));
sys = ss(Ac,B,C,0);
ic = [0,1,0];
t = linspace(0,3,100);
[x,t,y] = initial(sys,ic,t);
figure;
plot(t,x);

k1 = 50;
k2 = 10;
for i=1:n
    A = [0,1;-k1,-k2];
    B = [0,1]';
    C = [1,0];
    Q = eye(2);
    R = 1;
    [K,~,p] = lqr(A,B,Q,R);
    Ac = A-B*K;
    
    KK = -Ac(2,:);
    k1 = KK(1);
    k2 = KK(2);
    
end
disp(Ac);
disp(p);
fprintf('k1 = %f;\nk2 = %f;\n',KK(1),KK(2));
sys = ss(Ac,B,C,0);
ic = [1,0];
[x,t,y] = initial(sys,ic,t);
figure;
plot(t,x);