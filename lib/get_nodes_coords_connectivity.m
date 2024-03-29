function obj = get_nodes_coords_connectivity(dimensions, elements, dof_per_node)
obj.length = dimensions.length;
obj.width = dimensions.width;
obj.depth = dimensions.depth;
obj.total_elements = elements;
obj.total_nodes = elements + 1;
obj.total_dofs = obj.total_nodes * dof_per_node;
obj.connectivity = transpose([1:(obj.total_nodes-1);2:obj.total_nodes]);
obj.dof_per_element = dof_per_node * 2;
obj.global_coordinates = linspace(0,obj.length,obj.total_nodes);
obj.element_coordinates = zeros(obj.total_elements,2);
for i=1:obj.total_elements
    obj.element_coordinates(i,:) = obj.global_coordinates(obj.connectivity(i,:));
end
obj.dofs = node2dof(1:obj.total_nodes, dof_per_node);
if strcmp(dimensions.support_condition,'c')
    obj.arrested_nodes = 1;
    obj.arrested_dofs = node2dof(obj.arrested_nodes, dof_per_node);
    obj.arrested_dofs = obj.arrested_dofs(:);
    obj.free_dofs = setdiff(1:obj.total_dofs, obj.arrested_dofs);
elseif strcmp(dimensions.support_condition,'ss')
    obj.arrested_nodes = [1, obj.total_nodes];
    obj.arrested_dofs = node2dof(obj.arrested_nodes, dof_per_node);
    obj.arrested_dofs(end,:) = [];
    obj.arrested_dofs = obj.arrested_dofs(:);
    obj.free_dofs = setdiff(1:obj.total_dofs, obj.arrested_dofs);
end
end