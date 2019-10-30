close all;clear;clc;
addpath('lib');

dim.length = 1; % in m
dim.width = 0.1; % in m
dim.depth = 0.1; % in m
dim.support_condition = 'c'; % 'c' for cantilever

E_alum = 70e9; % in Pa
nu_alum = 0.35;
rho_alum = 2700; % in kg/m^3

aluminium = get_mechanical_properties(E_alum, nu_alum, rho_alum, -dim.width/2, dim.width/2, -dim.depth/2, dim.depth/2);

dof_per_node = 3;
elements = 10;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

force = [0;0;0];
point_load = -1;

max_iter = 20;
U = zeros(beam.total_dofs,1);
norm_residue = [];

udof = beam.dofs(1,:);
wdof = beam.dofs(2,:);
txdof = beam.dofs(3,:);

n_f_dof = beam.dofs(2,end);

tol_residue = 1e-8;

dU = zeros(beam.total_dofs, 1);
T_global = zeros(beam.total_dofs, beam.total_dofs);
R_global = zeros(beam.total_dofs,1);

for element = 1:beam.total_elements
    dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
%     el_connect = [dof_address(1,:),reshape(beam.dofs(2:end,1:2),[1,4])];
    el_connect = dof_address(:);
    
    el_u = U(dof_address(1,:));
    el_w = U(dof_address(2,:));
    el_tx = U(dof_address(3,:));
    
    element_matrices = get_element_matrices(aluminium, transpose(beam.element_coordinates(element,:)),el_u, el_w, el_tx, force);
    
    T_global(el_connect, el_connect) = T_global(el_connect, el_connect) + element_matrices.tangent;
    R_global(el_connect) = R_global(el_connect) + element_matrices.residue;
end
T = T_global(beam.free_dofs, beam.free_dofs);
R_global(n_f_dof) = point_load;
R = R_global(beam.free_dofs);

dU(beam.free_dofs) = T\R;
U = U + dU;