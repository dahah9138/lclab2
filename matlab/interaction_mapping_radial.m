clear
clf
clc

file = 'D:/lclab2 data/interactions/mat files/xhat_1const';
%[data, ~, ~] = load_interaction_data([file,'.bin']);
load([file,'.mat'],'data');

% Energy is expressed in units of kbT
energy = data(4,:);
r = data(3,:);

% Energy only matters to a factor of 1000 kbT
energy = round(energy);

scatter(r, (energy));
hold on
plot(r, (energy));
xlabel('Separation (p)')
ylabel('F (k_BT)')
axis([min(r) max(r) min(energy) max(energy)])