clear
clf
clc

file = 'D:/lclab2 data/interactions/mat files/sph_s2_5_3const';
half_sphere = 1;
upsample = 4;

load([file,'.mat'],'data');

% Throw out the first data point
data = data(:,2:end);

radius = 0.5*data(3,:);

% Energy is expressed in units of kbT
energy = round(data(4,:));

% F(phi, theta) = energy
F = scatteredInterpolant(data(2,:)', pi/2 - data(1,:)', energy');
pts = length(data(1,:));

[pq, tq] = meshgrid(linspace(0, 2*pi, pts*upsample),...
linspace(pi/2, 0, pts*upsample));


vq = F(pq, tq);

[X,Y,Z] = sph2cart(pq, tq, 1);

f1 = figure(1);
h = surf(X,Y,Z);
h.CData = vq;
set(h,'linestyle','none')
set(h, 'SpecularStrength', 0)
set(h, 'AmbientStrength', .5)

hold on

if half_sphere
    hbot = surf(-X,-Y,-Z);
    hbot.CData = vq;
    set(hbot,'linestyle','none')
    set(hbot, 'SpecularStrength', 0)
    set(hbot, 'AmbientStrength', .5)
end

axis equal off
grid off
%camlight

f1.Color = [1 1 1];
colorbar

f2 = figure(2);
% Polar sheet plot
Vq = interp2(vq,2);
[m,n] = size(Vq);

th = 180 / pi * linspace(pi/2,0,n);
ph = 180 / pi * linspace(0,2*pi,m);

imagesc(ph,th,Vq); colormap parula; axis xy;
colorbar;
xlabel('\phi (\circ)');
ylabel('\theta_{el} (\circ)');