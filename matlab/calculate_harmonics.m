clear
clc

file = 'D:/lclab2 data/heliknoton_new_3const';
R = 0.73; % Radius at which the harmonics are evaluated
component = 'x';
qn = 4; % The principle quantum number for which to compute harmonics
density = 50;
sampleThetaPoints = 25;
samplePhiPoints = 2*sampleThetaPoints;
thetaPoints = density;
phiPoints = 2*density;

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
v = zeros(sampleThetaPoints,samplePhiPoints);
thetaArr = linspace(0,pi,sampleThetaPoints);
phiArr = linspace(0,2*pi,samplePhiPoints);

% Assume the dimensions of the cube are [-1,1]^3
% Fill the function data to compute harmonics for
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
        v(t,p) = chi_tilt(xi,yi,zi);
    end
end

% Create the interpolant
[T,P] = ndgrid(thetaArr,phiArr);
F = griddedInterpolant(T,P,v);

% Create sphere to evaluate
evalThetaArr = linspace(0,pi,thetaPoints);
evalPhiArr = linspace(0,2*pi,phiPoints);
[evalT,evalP] = ndgrid(evalThetaArr,evalPhiArr);

% Compute spherical harmonic coefficients
coeffs = spherical_harmonics(evalT,evalP,F,qn);

% Print coefficients and find the dominant harmonics
l_dom = 0;
m_dom = 0;
ylm_dom_abs = 0;
for l=0:(qn-1)
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
% Evaluate the truncated series
[X,Y] = ndgrid(thetaArr,phiArr);
for i=1:qn
    res=f_Ylm(X,Y,coeffs,i);
    fprintf('Residual (n=%d): %f\n',i,sum(sum(abs(res-v))))
end

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
v = F(th,phi);
phase = (angle(v)/pi+1)/2;
C=hsv2rgb(phase,ones(size(phase)),ones(size(phase)));
[x,y,z] = sph2cart(phi,pi/2-th,abs(v));
s2=surf(x,y,z,abs(v));
s2.CData=C;
daspect([1 1 1])
