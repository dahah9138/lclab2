clear
clf
clc
fname = 'dz20';
floc = 'D:\dev\lclab2\data\interactions\';
[data, ~, ~] = load_interaction_data([floc,fname,'.bin']);
%load([file,'.mat'],'data');

% Energy is expressed in units of kbT
energy = data(4,:);
r = data(3,:);

% Energy only matters to a factor of 1 kbT
energy = round(energy);

f = figure(1);

scatter(r, (energy));
hold on
plot(r, (energy));
xlabel('Separation (p)')
ylabel('F (k_BT)')
xlim([min(r) max(r)])
ylim([min(energy) max(energy)])

%exportgraphics(f,[floc,fname,'.png'],'Resolution',300);