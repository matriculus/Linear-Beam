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
elements = 100;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

force = [0;0];
point_load = 0;

wdof = beam.dofs(1,:);
txdof = beam.dofs(2,:);

n_f_dof = beam.dofs(2,end);

element_act = [0];
element_sen = [0];

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
    l = beam.element_coordinates(element,:)*[-1;1];
    
    element_matrices = get_element_matrices(material.D, material.rho, l);
    
    K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
    M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
    F_global(el_connect,:) = F_global(el_connect,:) + element_matrices.force;
    
    K_beam(el_connect, el_connect) = K_beam(el_connect, el_connect) + element_matrices.stiffness;
    M_beam(el_connect, el_connect) = M_beam(el_connect, el_connect) + element_matrices.mass;
    
    for el = 1:n_act
        if any(element_act(el,:) == element)
            pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.piezoelectric_constant*beam.width, l, piezo_act.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            P_global(el_connect,el) = P_global(el_connect,el) + pzt_matrices.force;
        end
    end
    for el = 1:n_sen
        if any(element_sen(el,:) == element)
            pzt_matrices = get_sensor_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.piezoelectric_constant*pzt_depth/piezo_sen.dielectric_constant, l, piezo_sen.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            C_global(el,el_connect) = C_global(el,el_connect) + transpose(pzt_matrices.force);
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

%% Figures
shapes = zeros(beam.total_dofs,length(beam.free_dofs));
shapes(beam.free_dofs,:) = eigenvectors;

folder = sprintf('ModeFigures');
mkdir(folder);

for i=1:6
figure('position', [500, 500, 500, 100]);
w = shapes(wdof,i);
w = w/max(abs(w));
plot(beam.global_coordinates,w,'k','LineWidth',5);
% xlabel('X axis (m)');
% ylabel('W (m)');
axis([0,0.5,-1.2,1.2]);
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis
% tle = sprintf('Mode: %i',i);
% title(tle);

fname = sprintf('%s/Cantilever_beam_mode_%i',folder,i);
saveas(gcf,fname,'png');
saveas(gcf,fname,'fig');
end