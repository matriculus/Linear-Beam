close all;clear;clc;
addpath('lib');

%% Dimensions
dim.length = 1; % in m
dim.width = 0.01; % in m
dim.depth = 0.01; % in m
dim.support_condition = 'c'; % 'c' for cantilever

%% Material
E_alum = 70e9; % in Pa
nu_alum = 0.35;
rho_alum = 2700; % in kg/m^3

aluminium = get_mechanical_properties(E_alum, nu_alum, rho_alum, dim.width, -dim.depth/2, dim.depth/2);

%% Element properties
dof_per_node = 2;
elements = 10;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

force = [0;0;0];
point_load = -0.5;

wdof = beam.dofs(1,:);
txdof = beam.dofs(2,:);

n_f_dof = beam.dofs(2,end);

U = zeros(beam.total_dofs,1);
K_global = zeros(beam.total_dofs, beam.total_dofs);
M_global = zeros(beam.total_dofs, beam.total_dofs);
F_global = zeros(beam.total_dofs,1);

%% Matrix construction
for element = 1:beam.total_elements
    dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
    el_connect = dof_address(:);
    
    element_matrices = get_element_matrices(aluminium.D, aluminium.rho, beam.length);
    
    K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
    M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
    F_global(el_connect) = F_global(el_connect) + element_matrices.force;
end
K = K_global(beam.free_dofs, beam.free_dofs);
M = M_global(beam.free_dofs, beam.free_dofs);
F_global(n_f_dof) = point_load;
F = F_global(beam.free_dofs);

U(beam.free_dofs) = K\F;

n = length(beam.free_dofs);

Asys = [zeros(n,n), eye(n);-M\K, zeros(n,n)];
Bsys = [zeros(n,1); F];

tn = 200;
tspan = linspace(0,10,tn);
x0 = zeros(n, 1);
xd0 = zeros(n,1);
y0 = [x0;xd0];

ydot = @(t,y) Asys*y + Bsys*point_load;

[~, y] = ode45(ydot, tspan, y0);

U = zeros(tn, beam.total_dofs);
V = zeros(tn, beam.total_dofs);

U(:,beam.free_dofs) = y(:,1:n);
V(:,beam.free_dofs) = y(:,(n+1):2*n);

figure;
surf(U(:,wdof));