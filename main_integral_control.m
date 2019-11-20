close all;clear;clc;
addpath('lib');

%% Dimensions
dim.length = 0.5; % in m
dim.width = 0.03; % in m
dim.depth = 0.002; % in m
dim.support_condition = 'c'; % 'c' for cantilever

%% Material
% steel
E = 210e9; % in Pa
nu = 0.3;
rho = 7800; % in kg/m^3
cdr = 0.002; % critical damping ratio
material = get_mechanical_properties(E, nu, rho, dim.width, -dim.depth/2, dim.depth/2);

%% Actuator
% piezo
E_piezo = 139e9; % in Pa
nu_piezo = 0.3;
rho_piezo = 7500; % in kg/m^3

pzt_depth = 0.2e-3; % in m

piezo_act = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, dim.depth/2, dim.depth/2 + pzt_depth);
piezo_sen = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, -dim.depth/2-pzt_depth, -dim.depth/2);
piezo_act.piezoelectric_constant = 15.29; % Cm^-2
piezo_sen.piezoelectric_constant = 15.29; % Cm^-2
piezo_act.dielectric_constant = 11e-9; % Fm^-2
piezo_sen.dielectric_constant = 11e-9; % Fm^-2
%% Element
dof_per_node = 2;
elements = 10;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

force = [0;0];
point_load = 0;

wdof = beam.dofs(1,:);
txdof = beam.dofs(2,:);

n_f_dof = beam.dofs(2,end);

element_act = [1:3;4:6;7:9];
element_sen = [1:3;4:6;7:9];

n_act = size(element_act,1);
n_sen = size(element_sen,1);

%% Matrix construction
StaticU = zeros(beam.total_dofs,1);
K_global = zeros(beam.total_dofs, beam.total_dofs);
M_global = zeros(beam.total_dofs, beam.total_dofs);
F_global = zeros(beam.total_dofs,3);
P_global = zeros(beam.total_dofs,n_act);
C_global = zeros(n_sen,beam.total_dofs);

for element = 1:beam.total_elements
    dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
    el_connect = dof_address(:);
    x1 = beam.element_coordinates(element,1);
    x2 = beam.element_coordinates(element,2);
    
    element_matrices = get_element_matrices(material.D, material.rho, x1,x2);
    
    K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
    M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
    F_global(el_connect,:) = F_global(el_connect,:) + element_matrices.force;
    
    for el = 1:n_act
        if any(element_act(el,:) == element)
            pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.piezoelectric_constant,beam.width, x1, x2, piezo_act.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            P_global(el_connect,el) = P_global(el_connect,el) + pzt_matrices.force;
        end
    end
    for el = 1:n_sen
        if any(element_sen(el,:) == element)
            pzt_matrices = get_sensor_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.piezoelectric_constant,piezo_sen.dielectric_constant,pzt_depth, x1, x2, piezo_sen.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            C_global(el,el_connect) = C_global(el,el_connect) + pzt_matrices.force;
        end
    end
end
K = K_global(beam.free_dofs, beam.free_dofs);
M = M_global(beam.free_dofs, beam.free_dofs);
F = F_global(beam.free_dofs,:);
P = P_global(beam.free_dofs,:);
C = C_global(:,beam.free_dofs);

%% Static
StaticU(beam.free_dofs) = K\P*[1;-1;1];
figure;
plot(StaticU(wdof));

%% Modal analysis
[eigenvectors, eigenvalues] = eig(K,M);
[eigenvalues, indices] = sort(diag(eigenvalues));
eigenvectors = eigenvectors(:,indices);
frequencies = sqrt(eigenvalues);
w1 = frequencies(1);
w2 = frequencies(2);
bet = 2*cdr/(w1+w2);
alp = w1*w2*bet;

D = alp*M + bet*K;
%% Dynamics
n = length(beam.free_dofs);

Asys = [zeros(n,n), eye(n);-M\K, -M\D];
% Asys = [zeros(n,n), eye(n);-M\K, zeros(n,n)];
Bext = [zeros(n,3); M\F];
Bcont = [zeros(n,n_act); M\P];
Csys = [C, zeros(n_sen,n)];

tn = 100;
tf = 10;
tspan = linspace(0,tf,tn);
x0 = zeros(n,1);
xd0 = zeros(n,1);
y0 = [x0;xd0];

%% Controller
Q = 1*eye(2*n);
R = 0.01;
[Gain,~,~] = lqr(Asys,Bcont, Q, R);


%% Dynamic Inversion
shp = [1];
shape = zeros(n,1);
shape(beam.free_dofs) = eigenvectors(:,shp)*ones(length(shp),1);
shape = 1e-4*shape/max(shape(wdof));
xref = shape(beam.free_dofs);

kI = 200;
k1 = 100;
k2 = 10;
CMF = C*(M\P);
KG = [k1*C - C*(M\K), k2*C - C*(M\D), kI*eye(n_sen,n_sen)];
Gain1 = CMF\KG;
resF = CMF\(k1.*C)*xref;

%% Augmented system
y0 = [y0;zeros(n_sen,1)];
Aaug = [Asys,zeros(2*n,n_sen);Csys,zeros(n_sen,n_sen)];
Baug = [Bcont;zeros(n_sen,n_act)];
Be_aug = [Bext;zeros(n_sen,3)];
Ref = [zeros(2*n,n_sen);-eye(n_sen)]*C*xref;
%% ODE solution
f = @(t,ti) (t>ti).*(-1e-6.*sin(20.*(t-ti)) - 1e-6);
% ydot = @(t,y) Asys*y + Bcont*(-Gain)*y + Bext*[0;-1;0];
ydot = @(t,y) Aaug*y + Baug*(-Gain1)*y + Baug*resF + Be_aug*[0;f(t,0.0);0] + Ref;

[~, y] = ode45(ydot, tspan, y0);
%% Solution
U = zeros(tn, beam.total_dofs);
V = zeros(tn, beam.total_dofs);

U(:,beam.free_dofs) = y(:,1:n);
V(:,beam.free_dofs) = y(:,(n+1):2*n);

norm_con = zeros(tn,1);
for i=1:tn
    norm_con(i) = norm(U(i,:),2);
end
shape_norm = norm(shape,2).*ones(tn,1);
%% Input Output
Output = y(:,1:2*n)*Csys';
Input = y*-Gain1' + ones(tn,1)*resF';

%% Error in DI
Error = Output - ones(tn,1)*transpose(C*xref);
%% Figures
solution.Ucon = U;
solution.Vcon = V;
solution.Output = Output;
solution.Input = Input;
solution.norm_con = norm_con;
folder = sprintf('results/%s/forced_integral_control_[k1_%d,k2_%d,kI_%d]_s_%s_t_%0.1f',date,k1,k2,kI,num2str(shp),tf);
mkdir(folder);
fname = sprintf('%s/shape',folder);
filename = 'variables';
filename = sprintf('%s_%s.mat',fname,filename);
save(filename,'tspan','solution','beam');

lga = [];
lgs = [];
for i=1:n_act
    lga = [lga;sprintf('Actuator_%i',i)];
end

for i=1:n_sen
    lgs = [lgs;sprintf('Sensor_%i',i)];
end

figure;
surf(U(:,wdof));
xlabel('Length (m)');
ylabel('Time (s)');
zlabel('Deformation (m)');
filename = 'space_time';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

figure;
plot(tspan, Output,'LineWidth',2);
legend(lgs);
xlabel('Time (s)');
ylabel('Voltage (V)');
filename = 'sensor';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

figure;
plot(tspan, Input,'LineWidth',2);
legend(lga);
xlabel('Time (s)');
ylabel('Voltage (V)');
filename = 'actuator';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

figure;
plot(tspan, Error,'LineWidth',2);
xlabel('Time (s)');
ylabel('Error');
filename = 'error';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

figure;
plot(tspan, norm_con,'b-','LineWidth',2);
hold on;
plot(tspan, shape_norm,'r--','LineWidth',2);
hold off;
legend('Achieved','Target');
xlabel('Time (s)');
ylabel('Deformation norm (m)');
filename = 'norm_comp';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

%% shape comparison
figure;
plot(beam.global_coordinates,U(end,wdof)','b-','LineWidth',2);
hold on;
plot(beam.global_coordinates,shape(wdof),'r--','LineWidth',2);
legend('Achieved','Target');
xlabel('X axis (m)');
ylabel('Displacement (m)');
hold off;
filename = 'shape_comp';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');