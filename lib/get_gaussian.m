function obj = get_gaussian(coordinates, varargin)
I = ones(size(coordinates));
O = zeros(size(coordinates));
if isempty(varargin)
    obj.gp = int16(2);
else
    obj.gp = int16(varargin{1});
end
natural_coordinates = [-1;1];
if obj.gp==2
    Ig = ones(obj.gp,1);
    obj.points = [-0.5774; 0.5774];
    obj.weights = [1, 1];
    obj.matrix = [I,natural_coordinates];
    ply = [Ig, obj.points];
    obj.coordinates = (ply/obj.matrix)*coordinates;
end

end