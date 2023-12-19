clear
clf
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Files to load
fnames = ["field/interaction_field0_1","field/interaction_field0_4","field/interaction_field0_6"];
lgnd = ["E = 0.1", "E = 0.4", "E = 0.6"];
%% Parametrized path to trace
theta1 = 0;
dphi = 0*pi/180;
phi1 = pi/2-dphi;
theta2 = 30*pi/180;
phi2 = 3*pi/2+dphi;
tline = linspace(0,1,30);
param_path = zeros(30, 2);
% Straight line
%param_path(:,1) = (1-tline)*theta1 + tline*theta2;
%param_path(:,2) = (1-tline)*phi1 + tline*phi2;
% Parabolic path
param_path(:,2) = (1-tline)*phi1 + tline*phi2;
param_path(:,1) = theta2*(1-((param_path(:,2)-pi)).^2/(0.5*phi2-0.5*phi1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin the script
f1 = figure(1);
tc = linspace(0,1,length(lgnd));
% Choose colors with gray scale
c = ones(length(lgnd),3);
for i = 1:length(lgnd)
    c(i,:) = (1-tc(i)) * [.2,.2,.9] + tc(i) * [.9,.2,.2];
end
i = 1;
for fname=fnames
    F = load_interpolant(fname);
    param_field = F(param_path(:,2), param_path(:,1));
    % Shift to zero point
    param_field = param_field - param_field(1);
    plot(tline*180+90, param_field, 'LineWidth', 2,'Color', c(i,:));
    i = i + 1;
    hold on
end
legend(lgnd,'Location','north');
xlabel('\phi (\circ)');
ylabel('Free Energy (k_BT)');
xlim([90,270])
yticks(flip([0,-600,-1200,-1800]))
xticks([90,135,180,225,270])
hold off
exportgraphics(f1,strcat("D:/dev/lclab2/data/interactions/field/",strjoin(lgnd),"_multicut.png"),'Resolution',300)

function F = load_interpolant(fname)
    [data,~,~] = load_interaction_data(strcat("D:/dev/lclab2/data/interactions/",fname,".bin"));
    % Throw out the first data point
    data = data(:,2:end);
    % Energy is expressed in units of kbT
    energy = round(data(4,:));
    % F(phi, theta) = energy
    F = scatteredInterpolant(data(2,:)', pi/2 - data(1,:)', energy');
end

% Loads the new data format
function F = load_interpolant2(fname)
    [data,~,~] = load_interaction_data2(strcat("D:/dev/lclab2/data/interactions/",fname,".bin"));
    % Throw out the first data point
    data = data(:,2:end);
    % Energy is expressed in units of kbT
    energy = round(data(4,:));
    % F(phi, theta) = energy
    F = scatteredInterpolant(data(2,:)', pi/2 - data(1,:)', energy');
end
