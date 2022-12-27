clear
clc
addpath('./harmonic_lib/');
file = 'D:/lclab2 data/base_heliknoton';
R = 0.6; % Radius at which the harmonics are evaluated
component = 'x';
Nord = 3; % How many harmonics to compute
sampleThetaPoints = 60;
samplePhiPoints = sampleThetaPoints;
aziRes = 20;
polarRes = 20;
dirs = grid2dirs(aziRes, polarRes);

% Generate two random complex band-limited functions
Nf = 7;
Fnm = randn((Nf+1)^2,1) + 1i*randn((Nf+1)^2,1);
Ng = 5;
Gnm = randn((Ng+1)^2,1) + 1i*randn((Ng+1)^2,1);

% Generate a uniform grid for 12th-order integration
[~, tdesign_grid] = getTdesign(Nf+Ng);

% convert elevation to inclinations due to definition of spherical harmonics
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Read in the .lmt file
chi_file = [file,'_chi.mat'];
if isfile(chi_file)
    load(chi_file,'chi');
else
    data = fread_foel([file,'.lmt']);
    % Extract directors
    nn = data{10,2};
    % Compute the chi field
    chi = f_chi(nn, 0.05);
    save(chi_file,'chi');
end

% Make all chi points are tilted up
zhat = zeros(size(chi));
zhat(:,:,:,3) = 1;
chi(dot(chi,zhat,4) < 0) = -chi(dot(chi,zhat,4) < 0);

% Take x component of chi vector
if component == 'x'
    chi_tilt = chi(:,:,:,1);
else
    chi_tilt = chi(:,:,:,2);
end

[U,V,W] = size(chi_tilt);

% Find all points distance R from the center
thetaArr = linspace(0,pi,sampleThetaPoints);
phiArr = linspace(0,2*pi,samplePhiPoints);
regular_grid = grid2dirs(sampleThetaPoints,samplePhiPoints);
v = zeros(size(regular_grid(:,1)));
F_tdes = zeros(size(tdesign_grid(:,1)));
F = zeros(sampleThetaPoints,samplePhiPoints);
% Assume the dimensions of the cube are [-1,1]^3

for t=1:sampleThetaPoints
    for p=1:samplePhiPoints
        theta = thetaArr(t);
        phi = phiArr(p);

        % Get distance from center
        x = R * sin(theta) * cos(phi);
        y = R * sin(theta) * sin(phi);
        z = R * cos(theta);
        
        % Convert x,y,z to index space
        xi = int32((x+1)/2 * U);
        yi = int32((y+1)/2 * V);
        zi = int32((z+1)/2 * W);
        % Evaluate field data onto sphere
        F(t,p) = chi_tilt(xi,yi,zi);
    end
end

% Create the interpolant
[T,P] = ndgrid(thetaArr,phiArr);
F_interp = griddedInterpolant(T,P,F);

for t=1:length(regular_grid(:,1))
    theta = regular_grid(t,1);
    phi = regular_grid(t,2);

    % Evaluate field data onto sphere
    v(t) = F_interp(theta,phi);
end

for t=1:length(tdesign_grid(:,1))
    theta = tdesign_grid(t,1);
    phi = tdesign_grid(t,2);

    % Evaluate field data onto sphere
    F_tdes(t) = F_interp(theta,phi);
end

% get integration weights for the (non-uniform) regular grid, based on
% areas of voronoi cells around the sampling points
regular_weights = getVoronoiWeights(aziElev2aziIncl(regular_grid));
% perform SHT for the regular grid using weighted least-squares and complex SHs
coeffs = leastSquaresSHT(Nord, v, regular_grid, 'complex', regular_weights);
coeffs2 = directSHT(Nord, F_tdes, tdesign_grid, 'complex', []);

% Print coefficients and find the dominant harmonics
l_dom = 0;
m_dom = 0;
ylm_dom_abs = 0;
for l=0:Nord
    for m=-l:l
        idx = int32(l*(l+1) + m + 1);
        eval = abs(coeffs(idx));
        if eval > ylm_dom_abs
            l_dom = l;
            m_dom = m;
            ylm_dom_abs = eval;
        end
        fprintf('(l=%d, m=%d): |Y_lm| = %f\n',l,m,eval);
    end
end
fprintf('Dominant (l=%d, m=%d): |Y_lm| = %f\n',l_dom,m_dom,ylm_dom_abs);

% Show the dominant harmonic
th = linspace(0,pi,50);    % inclination
phi = linspace(0,2*pi,50); % azimuth
[th,phi] = meshgrid(th,phi);
% compute spherical harmonic
Y = harmonicY(l_dom,m_dom,th,phi,'type','real');
% plot the magnitude and magnitude
phase = (angle(Y)/pi+1)/2;
C=hsv2rgb(phase,ones(size(phase)),ones(size(phase)));
[x,y,z] = sph2cart(phi,pi/2-th,abs(Y));

figure(1)
s=surf(x,y,z,abs(Y));
s.CData=C;
daspect([1 1 1])

figure(2)
th = linspace(0,pi,50);    % inclination
phi = linspace(0,2*pi,50); % azimuth
[th,phi] = ndgrid(th,phi);
v = F_interp(th,phi);
phase = (angle(v)/pi+1)/2;
C=hsv2rgb(phase,ones(size(phase)),ones(size(phase)));
[x,y,z] = sph2cart(phi,pi/2-th,abs(v));
s2=surf(x,y,z,abs(v));
s2.CData=C;
daspect([1 1 1])

% recreate the function at a dense grid from the SH coefficients (using
% any of the above calculated coefficients) and plot reconstruction

figure(3)
plotSphFunctionCoeffs(coeffs, 'complex', 5, 5, 'real', gca); view(3), title('Reconstructed function')