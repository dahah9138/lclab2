clear
clf
clc

file = 'D:/lclab2 data/interactions/xhat_200_iter_20npp.bin';
[data, num_points, interaction_type] = load_interaction_data(file);

%file2 = 'D:/lclab2 data/interactions/zhat_200_iter_20npp_ext.bin';
%[data2, num_points, interaction_type] = load_interaction_data(file2);

%data = cat(2,data,data2);

energy = data(4,:);
en_sigma = data(5,:);
r = data(3,:);

% Choose the largest standard deviation
sigma = max(en_sigma);
en_precision = abs(floor(log10(sigma)));

% Energy only matters to a factor of 1000 kbT
energy = round(energy);

if isinf(en_precision) == false
    energy = round(energy, 4);
end

scatter(r, (energy));
hold on
plot(r, (energy));
xlabel('Separation (p)')
ylabel('F/V')