function dof = node2dof(nodes, dof_per_node)
n = length(nodes);
dof = zeros(dof_per_node, n);
for i = 1:dof_per_node
	dof(i,:) = nodes*dof_per_node - dof_per_node+i;
end
end