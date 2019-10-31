close all;clear;clc;
addpath('lib');

%% Dimensions
dim.length = 1; % in m
dim.width = 0.01; % in m
dim.depth = 0.01; % in m
dim.support_condition = 'c'; % 'c' for cantilever

%% Material properties
% aluminium
E_alum = 70e9; % in Pa
nu_alum = 0.35;
rho_alum = 2700; % in kg/m^3

aluminium = get_mechanical_properties(E_alum, nu_alum, rho_alum, dim.width, -dim.depth/2, dim.depth/2);


% piezo
E_piezo = 63e9; % in Pa
nu_piezo = 0.3;
rho_piezo = 7600; % in kg/m^3

depth = 0.2e-3; % in m

piezo_act = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, dim.depth/2, dim.depth/2 + depth);
piezo_sen = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, -dim.depth/2-depth, -dim.depth/2);
piezo_act.piezoelectric_constant = 17.584; % Cm^-2
piezo_sen.piezoelectric_constant = 17.584; % Cm^-2
piezo_act.dielectric_constant = 15e-9; % Fm^-2
%% Element properties
dof_per_node = 3;
elements = 20;
beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);

element_act = [1:elements];
element_sen = [1];

%% Load
force = [0;0;0];
point_load = 0;
voltage = 1000;
%% Initialisation
max_iter = 20;
U = zeros(beam.total_dofs,1);
norm_residue = [];

udof = beam.dofs(1,:);
wdof = beam.dofs(2,:);
txdof = beam.dofs(3,:);

n_f_dof = beam.dofs(2,end);

tol_residue = 1e-2;

q = zeros(max_iter,1);
ww = zeros(max_iter, 1);

for iter = 1:max_iter
    P = force * iter/max_iter;
    pt = point_load*iter/max_iter;
    V = voltage*iter/max_iter;
    q(iter) = V;
    residue_error = 1;
    j=0;
    %% Convergence
    while (residue_error > tol_residue && j < 25)
        %     for i=1:5
        j = j + 1;
        dU = zeros(beam.total_dofs, 1);
        T_global = zeros(beam.total_dofs, beam.total_dofs);
        R_global = zeros(beam.total_dofs,1);
        F_global = zeros(beam.total_dofs,1);
        
        for element = 1:beam.total_elements
            dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
            %             el_connect = [dof_address(1,:),reshape(beam.dofs(2:end,1:2),[1,4])];
            el_connect = dof_address(:);
            
            el_u = U(dof_address(1,:));
            el_w = U(dof_address(2,:));
            el_tx = U(dof_address(3,:));
            
            element_matrices = get_element_matrices(aluminium, transpose(beam.element_coordinates(element,:)),el_u, el_w, el_tx, P);
            
            T_global(el_connect, el_connect) = T_global(el_connect, el_connect) + element_matrices.tangent;
            R_global(el_connect) = R_global(el_connect) + element_matrices.residue;
            if any(element_act==element)
                pzt_matrices = get_actuator_matrices(piezo_act, transpose(beam.element_coordinates(element,:)),el_u, el_w, el_tx, V);
                T_global(el_connect, el_connect) = T_global(el_connect, el_connect) + pzt_matrices.tangent;
                F_global(el_connect) = F_global(el_connect) + pzt_matrices.residue;
            end
        end
        T = T_global(beam.free_dofs, beam.free_dofs);
        R = R_global(beam.free_dofs) + F_global(beam.free_dofs);
        
        dU(beam.free_dofs) = T\R;
        U = U + dU;
        residue_error = norm(dU)/norm(U);
        norm_residue = [norm_residue, residue_error];
        
        fprintf('Residue: %d\tRank(T): %i\n',residue_error, rank(T));
    end
    %%
    ww(iter) = max(abs(U(wdof)));
    fprintf('Iteration: %i\tDisplacement: %f mm\n',iter,ww(iter)*1e3);
    plot(U(wdof));
    hold on;
end

hold off;
