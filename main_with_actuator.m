close all;clear;clc;
addpath('lib');

%% Dimensions
dim.length = 1; % in m
dim.width = 0.1; % in m
dim.depth = 0.1; % in m
dim.support_condition = 'c'; % 'c' for cantilever

%% Material
E_alum = 70e9; % in Pa
nu_alum = 0.35;
rho_alum = 2700; % in kg/m^3

aluminium = get_mechanical_properties(E_alum, nu_alum, rho_alum, dim.width, -dim.depth/2, dim.depth/2);

%% Actuator
E_piezo = 63e9; % in Pa
nu_piezo = 0.3;
rho_piezo = 7600; % in kg/m^3

depth = 0.2e-3; % in m

piezo_act = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, dim.depth/2, dim.depth/2 + depth);
piezo_sen = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, -dim.depth/2-depth, -dim.depth/2);
piezo_act.e31 = 17.584; % Cm^-2
piezo_sen.e31 = 17.584; % Cm^-2
piezo_act.dielectric_constant = 15e-9; % Fm^-2

%% Element
dof_per_node = 2;
elements = 10;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

force = [0;0;0];
point_load = 0;

wdof = beam.dofs(1,:);
txdof = beam.dofs(2,:);

n_f_dof = beam.dofs(2,end);

element_act = [1:3];
element_sen = [0];

%% Matrix construction
U = zeros(beam.total_dofs,1);
K_global = zeros(beam.total_dofs, beam.total_dofs);
M_global = zeros(beam.total_dofs, beam.total_dofs);
F_global = zeros(beam.total_dofs,1);
C_global = zeros(beam.total_dofs,1);

for element = 1:beam.total_elements
    dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
    el_connect = dof_address(:);
    
    element_matrices = get_element_matrices(aluminium.D, aluminium.rho, beam.length);
    
    K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
    M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
    F_global(el_connect) = F_global(el_connect) + element_matrices.force;
    
    if any(element_act == element)
        pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.e31*beam.width, beam.length, piezo_act.lever_arm);
        K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
        M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
        F_global(el_connect) = F_global(el_connect) + pzt_matrices.force;
    end
    if any(element_sen == element)
        pzt_matrices = get_actuator_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.e31*beam.width, beam.length, piezo_sen.lever_arm);
        K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
        M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
        C_global(el_connect) = C_global(el_connect) + pzt_matrices.force;
    end
end
K = K_global(beam.free_dofs, beam.free_dofs);
M = M_global(beam.free_dofs, beam.free_dofs);
F = F_global(beam.free_dofs);

%% Static
U(beam.free_dofs) = K\F;

%% Modal analysis
[eigenvectors, eigenvalues] = eig(K,M);
[eigenvalues, indices] = sort(diag(eigenvalues));
eigenvectors = eigenvectors(:,indices);
frequencies = sqrt(eigenvalues);
%% Dynamics
n = length(beam.free_dofs);

Asys = [zeros(n,n), eye(n);-M\K, zeros(n,n)];
Bsys = [zeros(n,1); F];

tn = 200;
tspan = linspace(0,10,tn);
x0 = zeros(n, 1);
xd0 = zeros(n,1);
y0 = [x0;xd0];

ydot = @(t,y) Asys*y + Bsys*1;

[~, y] = ode45(ydot, tspan, y0);

U = zeros(tn, beam.total_dofs);
V = zeros(tn, beam.total_dofs);

U(:,beam.free_dofs) = y(:,1:n);
V(:,beam.free_dofs) = y(:,(n+1):2*n);

figure;
surf(U(:,wdof));