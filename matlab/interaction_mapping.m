clear
clf
clc

file = 'D:/lclab2 data/interactions/sph_s2_2_200_iter_20npp.bin';
half_sphere = 1;
upsample = 4;
% don't touch
s = dir(file);
double_sz = 8;
datastruct_size = 5 * double_sz;
num_points = (s.bytes-1) / datastruct_size;

fID = fopen(file, 'r');
interaction_type = fread(fID, 1, 'char*1');
data = fread(fID, [5 num_points], 'double');
fclose(fID);

% Throw out the first data point
data = data(:,2:end);
num_points = num_points - 1;

scan = [];

for i=1:num_points
    if abs(data(3,i)) < 1e6
        scan(:,end+1) = data(:,i);
    end
end

data = scan;
clear scan

radius = 0.5*data(3,:);

% Choose the largest standard deviation
energy = round(data(4,:));
en_sigma = data(5,:);
sigma = max(en_sigma);
en_precision = abs(floor(log10(sigma)));
if isinf(en_precision) == false
    %energy = round(energy, en_precision);
end

% F(phi, theta) = energy
F = scatteredInterpolant(data(2,:)', pi/2 - data(1,:)', energy');
pts = num_points;

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

fig = gcf;
fig.Color = [1 1 1];
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