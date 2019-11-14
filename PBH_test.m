% close all;
% %%
% close all;clear;clc;
% addpath('lib');
% 
% %% Dimensions
% dim.length = 1; % in m
% dim.width = 0.1; % in m
% dim.depth = 0.01; % in m
% dim.support_condition = 'c'; % 'c' for cantilever
% 
% %% Material
% E_alum = 70e9; % in Pa
% nu_alum = 0.35;
% rho_alum = 2700; % in kg/m^3
% cdr = 0.002; % critical damping ratio
% aluminium = get_mechanical_properties(E_alum, nu_alum, rho_alum, dim.width, -dim.depth/2, dim.depth/2);
% 
% %% Actuator
% E_piezo = 63e9; % in Pa
% nu_piezo = 0.3;
% rho_piezo = 7600; % in kg/m^3
% 
% depth = 1e-3; % in m
% 
% piezo_act = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, dim.depth/2, dim.depth/2 + depth);
% piezo_sen = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, -dim.depth/2-depth, -dim.depth/2);
% piezo_act.e31 = 17.584; % Cm^-2
% piezo_sen.e31 = 17.584; % Cm^-2
% piezo_act.dielectric_constant = 15e-9; % Fm^-2
% 
% %% Element
% dof_per_node = 2;
% elements = 10;
% beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);
% 
% force = [0;0];
% point_load = 0;
% 
% wdof = beam.dofs(1,:);
% txdof = beam.dofs(2,:);
% 
% n_f_dof = beam.dofs(2,end);
% 
% % test_a = [1:elements]';
% Sa_diag = [];
% % sa_small = zeros(size(test_a,1),1);
% % for ppz=1:size(test_a,1)
% % element_act = test_a(ppz,:);
% % element_sen = [0];
% 
% element_act = [1:3];
% element_sen = [1:3];
% 
% n_act = size(element_act,1);
% n_sen = size(element_sen,1);
% 
% %% Matrix construction
% 
% K_global = zeros(beam.total_dofs, beam.total_dofs);
% M_global = zeros(beam.total_dofs, beam.total_dofs);
% F_global = zeros(beam.total_dofs,3);
% P_global = zeros(beam.total_dofs,n_act);
% C_global = zeros(n_sen,beam.total_dofs);
% 
% for element = 1:beam.total_elements
%     dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
%     el_connect = dof_address(:);
%     l = beam.element_coordinates(element,:)*[-1;1];
%     
%     element_matrices = get_element_matrices(aluminium.D, aluminium.rho, l);
%     
%     K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
%     M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
%     F_global(el_connect,:) = F_global(el_connect,:) + element_matrices.force;
%     
%     for el = 1:n_act
%         if any(element_act(el,:) == element)
%             pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.e31*beam.width, l, piezo_act.lever_arm);
%             K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
%             M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
%             P_global(el_connect,el) = P_global(el_connect,el) + pzt_matrices.force;
%         end
%     end
%     for el = 1:n_sen
%         if any(element_sen(el,:) == element)
%             pzt_matrices = get_actuator_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.e31*beam.width, l, piezo_sen.lever_arm);
%             K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
%             M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
%             C_global(el,el_connect) = C_global(el,el_connect) + transpose(pzt_matrices.force);
%         end
%     end
% end
% K = K_global(beam.free_dofs, beam.free_dofs);
% M = M_global(beam.free_dofs, beam.free_dofs);
% F = F_global(beam.free_dofs,:);
% P = P_global(beam.free_dofs,:);
% C = C_global(:,beam.free_dofs);
% 
% %% Modal analysis
% [eigenvectors, eigenvalues] = eig(K,M);
% [eigenvalues, indices] = sort(diag(eigenvalues));
% eigenvectors = eigenvectors(:,indices);
% frequencies = sqrt(eigenvalues);
% w1 = frequencies(1);
% w2 = frequencies(2);
% bet = 2*cdr/(w1+w2);
% alp = w1*w2*bet;
% 
% D = alp*M + bet*K;
% %% Dynamics
% n = length(beam.free_dofs);
% 
% Asys = [zeros(n,n), eye(n);-M\K, -M\D];
% % Asys = [zeros(n,n), eye(n);-M\K, zeros(n,n)];
% Bext = [zeros(n,3); M\F];
% Bcont = [zeros(n,n_act); M\P];
% Csys = [C, zeros(n_sen,n)];
%%
% lmda = eig(Asys);
% rnk_pbh = zeros(size(lmda));
% for i = 1:length(lmda)
%     rnk_pbh(i) = rank([Asys - lmda(i)*eye(2*n), Bcont],1e-5);
% end
% 
% [Va, Da] = eig(Asys);
% 
% %% Gramian
% sys = ss(Asys,Bcont,Csys,0);
% WC = gram(sys,'c');
% WO = gram(sys,'o');
% [Ua,Sa,~] = svd(WC);
% Sa_diag = diag(Sa);
% % sa_small(ppz) = min(abs(diag(Sa)));
% % % end
% % plot(sa_small);
% qq = diag(y_c*WC*y_c');
% plot(qq)

%%
qq = rere(1);
function q = rere(w, varargin)

q = w;
end