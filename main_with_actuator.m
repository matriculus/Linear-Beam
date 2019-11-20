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
element_sen = element_act;

n_act = size(element_act,1);
n_sen = size(element_sen,1);

%% Matrix construction
K_global = zeros(beam.total_dofs, beam.total_dofs);
M_global = zeros(beam.total_dofs, beam.total_dofs);
F_global = zeros(beam.total_dofs,3);
P_global = zeros(beam.total_dofs,n_act);
C_global = zeros(n_sen,beam.total_dofs);
K_beam = zeros(beam.total_dofs, beam.total_dofs);
M_beam = zeros(beam.total_dofs, beam.total_dofs);

for element = 1:beam.total_elements
    dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
    el_connect = dof_address(:);
    x1 = beam.element_coordinates(element,1);
    x2 = beam.element_coordinates(element,2);
    
    element_matrices = get_element_matrices(material.D, material.rho, x1,x2);
    
    K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
    M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
    F_global(el_connect,:) = F_global(el_connect,:) + element_matrices.force;
    
    K_beam(el_connect, el_connect) = K_beam(el_connect, el_connect) + element_matrices.stiffness;
    M_beam(el_connect, el_connect) = M_beam(el_connect, el_connect) + element_matrices.mass;
    
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

Kb = K_beam(beam.free_dofs, beam.free_dofs);
Mb = M_beam(beam.free_dofs, beam.free_dofs);

%% Modal analysis
[eigenvectors, eigenvalues] = eig(Kb,Mb);
[eigenvalues, indices] = sort(diag(eigenvalues));
eigenvectors = eigenvectors(:,indices);
frequencies = sqrt(eigenvalues);
w1 = frequencies(1)/(2*pi);
w2 = frequencies(2)/(2*pi);
bet = 2*cdr/(w1+w2);
alp = w1*w2*bet;

D = alp*M + bet*K;

%% Frequencies
fprintf('Frequencies:\n');
for f=frequencies
    fprintf('%f\n',f);
end
%% Dynamics
n = length(beam.free_dofs);

Asys = [zeros(n,n), eye(n);-M\K, -M\D];
% Asys = [zeros(n,n), eye(n);-M\K, zeros(n,n)];
Bext = [zeros(n,3); M\F];
Bcont = [zeros(n,n_act); M\P];
Csys = [C, zeros(n_sen,n)];

tn = 100;
tf = 4;
tspan = linspace(0,tf,tn);
x0 = zeros(n,1);
xd0 = zeros(n,1);
y0 = [x0;xd0];

%% Controller
Q = 1*eye(2*n);
R = 0.01;
[Gain,~,~] = lqr(Asys,Bcont, Q, R);


%% Dynamic Inversion
shp = [2];
shape = zeros(n,1);
shape(beam.free_dofs) = eigenvectors(:,shp)*ones(length(shp),1);
shape = 1e-4*shape/max(shape(wdof));
xref = shape(beam.free_dofs);

k1 = 500;
k2 = 10;
CMF = C*(M\P);
KG = [k1*C - C*(M\K), k2*C - C*(M\D)];
Gain1 = CMF\KG;
resF = CMF\(k1.*C)*xref;

%% ODE solution
f = @(t,ti) (t>ti).*(-1e-6.*sin(20.*(t-ti)) - 1e-6);
ydot_c = @(t,y) Asys*y + Bcont*(-Gain1*y + resF) + Bext*[0;0.*f(t,0.0);0]; % (t>0.5)*1e-4*exp(-10*(t-0.5))

[~, y_c] = ode45(ydot_c, tspan, y0);

%% Post processing
Ucon = zeros(tn, beam.total_dofs);
Vcon = zeros(tn, beam.total_dofs);

Ucon(:,beam.free_dofs) = y_c(:,1:n);
Vcon(:,beam.free_dofs) = y_c(:,(n+1):2*n);

norm_con = zeros(tn,1);
for i=1:tn
    norm_con(i) = norm(Ucon(i,:),2);
end
shape_norm = norm(shape,2).*ones(tn,1);
%% Input Output
Output = y_c*Csys';
Input = y_c*-Gain1' + ones(tn,1)*resF';

%% Error in DI
Error = Output - ones(tn,1)*transpose(C*xref);
%% Figures
solution.Ucon = Ucon;
solution.Vcon = Vcon;
solution.Output = Output;
solution.Input = Input;
solution.norm_con = norm_con;
%%
folder = sprintf('results/%s/forced_vib_control_[k1_%d,k2_%d]_s_%s_t_%0.1f',date,k1,k2,num2str(shp),tf);
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
surf(real(Ucon(:,wdof)));
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
filename = 'Error';
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

figure;
plot(beam.global_coordinates, shape(wdof),'r--','LineWidth',2);
hold on;
plot(beam.global_coordinates, Ucon(end,wdof),'b-','LineWidth',2);
hold off;
legend('Target','Achieved');
xlabel('Length (m)');
ylabel('Displacement (m)');
filename = 'shape';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');
