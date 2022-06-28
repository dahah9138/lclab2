file = '../data/mat/Q1.mat';
load(file, 'nn');

[m,n,p,~] = size(nn);

% Top surface topography
z_profile = zeros(m, n);
h = zeros(m, n);
phi = zeros(m, n);

h_par = 0.9;
h0 = 1;
h_perp = 1.1;

for i=1:p
    % Determine the contribution to z_profile by nz
    phi = acos(nn(:,:,i, 3));
    h = sqrt((h_par*cos(phi).^2+(h_perp*sin(phi)).^2));
    z_profile = z_profile + h;
end

f = figure(1);

% Draw the z profile
s = surf(z_profile);
s.EdgeColor = 'none';

set(gcf,'Color','w')
axis vis3d
axis off
colorbar