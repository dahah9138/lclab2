% Plot the dominant harmonics for heliknoton attractions
th = linspace(0,pi,50);    % inclination
phi = linspace(0,2*pi,50); % azimuth
[th,phi] = meshgrid(th,phi);
% compute spherical harmonic
Y = harmonicY(2,1,th,phi,'type','real');
% plot the magnitude and magnitude
phase = (angle(Y)/pi+1)/2;
C=hsv2rgb(phase,ones(size(phase)),ones(size(phase)));
[x,y,z] = sph2cart(phi,pi/2-th,abs(Y));
s=surf(x,y,z,abs(Y));
s.CData=C;
hold on
Y = harmonicY(2,-1,th,phi,'type','real');
% plot the magnitude and magnitude
phase = (angle(Y)/pi+1)/2;
C=hsv2rgb(phase,ones(size(phase)),ones(size(phase)));
[x,y,z] = sph2cart(phi,pi/2-th,abs(Y));
s2=surf(x,y,z,abs(Y));
s2.CData=C;

daspect([1 1 1])