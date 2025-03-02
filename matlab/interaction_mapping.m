clear
clf
clc

%file = 'D:/lclab2 data/interactions/mat files/sph_s2_5_3const';
half_sphere = 1;
upsample = 1;
%fp = "field_final/";
%fname = "E0p6";
fname = "KnotInteraction_2025";
draw_arrows = 0;
exp_traj = 0;
use_contours = 0;
parabola_cut = 0;
%fname = "field/interaction_field0_1";
%load([file,'.mat'],'data');
% Interaction energy at zero applied voltage
%[data,~,~] = load_interaction_data('D:\dev\lclab2\data\interactions\interaction_E0.bin');
[data,~,~] = load_interaction_data2(strcat("D:/dev/lclab2/data/interactions/",fname,".bin"));

% Throw out the first data point
d_original = data;
d0 = data(:,1);
data = data(:,2:end);

radius = 0.5*data(3,:);

if draw_arrows == 1
    rvec1 = data(5:7,:);
    rvec2 = data(8:10,:);

    % Get the original vecs
    rhat(1,:) = sin(data(1,:)).*cos(data(2,:));
    rhat(2,:) = sin(data(1,:)).*sin(data(2,:));
    rhat(3,:) = cos(data(1,:));
    rvec1_i(1,:) = 0.5*data(3,:).*rhat(1,:);
    rvec1_i(2,:) = 0.5*data(3,:).*rhat(2,:);
    rvec1_i(3,:) = 0.5*data(3,:).*rhat(3,:);
    rvec2_i(1,:) = -rvec1_i(1,:);
    rvec2_i(2,:) = -rvec1_i(2,:);
    rvec2_i(3,:) = -rvec1_i(3,:);
    
    dr1 = rvec1 - rvec1_i;
    dr2 = rvec2 - rvec2_i;
    
    % Subtract off the component parallel to the radial vector
    dr1 = dr1 - dot(dr1,rhat,1).*rhat;
    dr2 = dr2 - dot(dr2,rvec2_i,1).*rvec2_i;

    % Normalization by the initial separation
    dr1 = dr1./(2*vecnorm(rvec1_i,2,1));
    dr2 = dr2./(2*vecnorm(rvec1_i,2,1));
end


% Energy is expressed in units of kbT
energy = round(data(4,:));

% F(phi, theta) = energy
F = scatteredInterpolant([d0(2),data(2,:)]', pi/2 - [d0(1),data(1,:)]', [d0(4),energy]');
pts = length(data(1,:));
dim = round(0.5+0.5*sqrt(1+4*pts));

[pq, tq] = meshgrid(linspace(0, 2*pi, dim*upsample),...
linspace(pi/2, 0, dim*upsample));
vq = F(pq, tq);

%% Original data map
% (phi1,theta1,E11), (phi2,theta1,E21), ... (phiN,theta1,EN1)
% (phi1,theta2,E12),...
% ...                                       (phiN,thetaN,ENN)
pmap = linspace(0,2*pi,dim);
tmap = linspace(0,pi/2,dim);
elemap = pi/2 - tmap;

Emap = zeros(dim,dim);
for t = 1:dim
    for p = 1:dim
        id = p + (t - 2) * dim;
        id
        if t == 1 % Degenerate point
            Emap(t,p) = d0(4);
        else % Nondegenerate points
            Emap(t,p) = data(4,id);
        end
    end
end

[pg,tg] = meshgrid(pmap,elemap);

[X,Y,Z] = sph2cart(pg, tg, 1);

f1 = figure(1);
h = surf(X,Y,Z);
h.CData = Emap;
set(h,'linestyle','none')
set(h, 'SpecularStrength', 0)
set(h, 'AmbientStrength', .5)

hold on

if half_sphere
    hbot = surf(-X,-Y,-Z);
    hbot.CData = Emap;
    set(hbot,'linestyle','none')
    set(hbot, 'SpecularStrength', 0)
    set(hbot, 'AmbientStrength', .5)
end

if draw_arrows == 1
    % Plots the displacement on the sphere surface with red arrows
    % (perpendicular components has been subtracted off)
    r0 = 1.02;
    hold on
    quiver3(r0*rhat(1,:),r0*rhat(2,:),r0*rhat(3,:),dr1(1,:),dr1(2,:),dr1(3,:),'r')
    hold on
    quiver3(-r0*rhat(1,:),-r0*rhat(2,:),-r0*rhat(3,:),dr2(1,:),dr2(2,:),dr2(3,:),'r')
end

axis equal off
grid off
%camlight

f1.Color = [1 1 1];
colorbar

f2 = figure(2);
% Polar sheet plot
Vq = interp2(Emap,3);
[m,n] = size(Vq);

ele = 180 / pi * linspace(pi/2,0,n);
ph = 180 / pi * linspace(0,2*pi,m);

imagesc(ph,ele,Vq); colormap parula; axis xy;
colorbar;
xlabel('\phi (\circ)');
ylabel('\theta_{el} (\circ)');

if use_contours
    hold on
    [mph,mth] = meshgrid(ph,ele);
    [~,cntour] = contour(mph,mth,Vq,[600,300,0,-300],'LineColor','k');
    %cntour.EdgeColor = [0 0 0];
    cntour.LineWidth = 3;
end
xlabel('\phi (\circ)','FontSize',20);
ylabel('\theta_{el} (\circ)','FontSize',20);

ax = gca;
ax.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 20;
ax.YLabel.FontSize = 25;

if exp_traj == 1
    hold on
    exp_trajectory('D:\Smalyukh\Hybridized heliknoton publication\Figures\figure_2\parameters.xlsx');
end

% Parametrize a linear line from theta = 0, phi = 90 deg to
% theta = 40 deg, phi = 180
theta1 = 0;
dphi = 0*pi/180;
phi1 = pi/2-dphi;
theta2 = 32*pi/180;
phi2 = 3*pi/2+dphi;
tline = linspace(0,1,30);
param_path = zeros(30, 2);

% Straight line
%param_path(:,1) = (1-tline)*theta1 + tline*theta2;
%param_path(:,2) = (1-tline)*phi1 + tline*phi2;
% Parabolic path
param_path(:,2) = (1-tline)*phi1 + tline*phi2;
param_path(:,1) = theta2*(1-((param_path(:,2)-pi)).^2/(0.5*phi2-0.5*phi1)^2);

if parabola_cut == 1
hold on
plot(param_path(:,2)*180/pi,param_path(:,1)*180/pi, 'r', 'LineStyle','--', 'LineWidth', 3);
end

param_field = F(param_path(:,2), param_path(:,1));

%f3 = figure(3);
%plot(tline, param_field,'k', 'LineWidth', 2);
%mkdir(strcat("D:/dev/lclab2/data/interactions/",fp,fname,"/"));
%exportgraphics(f1,strcat("D:/dev/lclab2/data/interactions/",fp,fname,"/",fname,"_sphere.png"),'Resolution',300)
%exportgraphics(f2,strcat("D:/dev/lclab2/data/interactions/",fp,fname,"/",fname,"_sheet.png"),'Resolution',300)
%exportgraphics(f3,strcat("D:/dev/lclab2/data/interactions/",fp,fname,"/",fname,"_cut.png"),'Resolution',300)