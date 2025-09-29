clear
clf
clc

%file = 'D:/lclab2 data/interactions/mat files/sph_s2_5_3const';
upsample = 4;
%fp = "field_final/";
%fname = "E0p6";
fname = "KnotInteraction_2025";
%load([file,'.mat'],'data');
% Interaction energy at zero applied voltage
%[data,~,~] = load_interaction_data('D:\dev\lclab2\data\interactions\interaction_E0.bin');
[data,~,~] = load_interaction_data2(strcat("D:/dev/lclab2/data/interactions/",fname,".bin"));

% Add degenerate points for theta = 0
data_original = data;
nTheta = 20;
nPhi = 40;
data_new = zeros(4,nTheta*nPhi);
data_new(:,1:nPhi) = repmat(data_original(:,1),1,nPhi);
data_new(:,(nPhi+1):end) = data_original(:,2:end);
data = data_new;

%data = data(:,2:end);
radius = 0.5*data(3,:);
theta_1d = pi/2-data(1,:);
phi_1d = data(2,:);
% Energy is expressed in units of kbT (T = 293 K)
energy_1d = round(data(4,:));

% Transform 1d data to 2d
theta_2d = reshape(theta_1d,[nPhi,nTheta])';
phi_2d = reshape(phi_1d,[nPhi,nTheta])';
phi_2d(1,:) = phi_2d(2,:);
energy_2d = reshape(energy_1d,[nPhi,nTheta])';

ph = 180/pi*[phi_2d(1,1),phi_2d(1,end)];
th = 180/pi*[theta_2d(1,1),theta_2d(end,1)];

% Interp the grid if upsampling
if upsample > 1
    energy_2d = interp2(energy_2d,upsample);
end


imagesc(ph,th,energy_2d); colormap parula; axis xy;
colorbar;
ax = gca;
ax.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 20;
ax.YLabel.FontSize = 25;
xticks([0,90,180,270,360]);
xlabel('$\phi$ ($^\circ$)','FontSize',20,'Interpreter','latex');
ylabel('$\theta_{el}$ ($^\circ$)','FontSize',20,'Interpreter','latex');
