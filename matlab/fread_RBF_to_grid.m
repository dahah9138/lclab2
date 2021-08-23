function [nn, vv] = fread_RBF_to_grid(file)

data = fread_RBF(file);
cdims = cell2mat(data(2,2));
directors = cell2mat(data(3,2));
positions = cell2mat(data(4,2));
voltage = cell2mat(data(5,2));
active_nodes = cell2mat(data(6,2));

n_active_scattered = directors(active_nodes,:);
p_active_scattered = positions(active_nodes,:);
v_active_scattered = voltage(active_nodes,:);

% Three interpolant functions to cast to grid

Fx = scatteredInterpolant(p_active_scattered, n_active_scattered(:,1));
Fy = scatteredInterpolant(p_active_scattered, n_active_scattered(:,2));
Fz = scatteredInterpolant(p_active_scattered, n_active_scattered(:,3));

V = scatteredInterpolant(p_active_scattered, v_active_scattered);

npp = 20;
cdimshalf = cdims ./ 2;
m = cdims(1)*npp;
n = cdims(2)*npp;
p = cdims(3)*npp;
x = linspace(-cdimshalf(1),cdimshalf(1), m);
y = linspace(-cdimshalf(2),cdimshalf(2), n);
z = linspace(-cdimshalf(3),cdimshalf(3), p);

% Assemble 3D grid
[xx, yy, zz] = ndgrid(x,y,z);

% Sample F at the points

disp("Applying interpolant")

nx = Fx(xx,yy,zz);

disp("nx complete")

ny = Fy(xx,yy,zz);

disp("ny complete")

nz = Fz(xx,yy,zz);

disp("nz complete")

nn = cat(4, nx, ny, nz);

vv = V(xx,yy,zz);

disp("v complete")

end