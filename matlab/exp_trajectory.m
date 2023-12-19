function exp_trajectory(file)

data = readtable(file);

x1 = data.x5;
y1 = data.y5;
x2 = data.x6;
y2 = data.y6;

x = x2 - x1;
y = y2 - y1;
p = 71.5/2; % Pitch in units pixels (half heliknoton long axis)

psi = data.psi56/180*pi;
phi = atan2(y,x) + pi/2;

rho = sqrt(x.*x + y.*y);
% Throw away points with rho = 0
psi = psi(rho > 0);
phi = phi(rho > 0);
rho = rho(rho > 0);

theta_el = asin(psi./(2*pi)*p./sqrt(rho.^2+(p*psi/2/pi).^2));

% Plot the trajectory taken
plot(180/pi*phi,180/pi*theta_el, 'r', 'LineWidth',3)

end