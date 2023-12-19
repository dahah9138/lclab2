clear
clf
clc

dz_array = [6,7,8,9,10,15,20];
amplitude = zeros(size(dz_array));

for i=1:length(dz_array)

    floc = "D:\dev\lclab2\data\interactions\dz";
    fname = strcat(floc,num2str(dz_array(i)),".bin")
    [data, ~, ~] = load_interaction_data(fname);
    %load([file,'.mat'],'data');
    
    % Energy is expressed in units of kbT
    energy = data(4,:);
    r = data(3,:);
    
    % Energy only matters to magnitude of 1 kbT
    energy = round(energy);

    amplitude(i) = max(energy) - min(energy);

end

f = figure(1);

scatter(dz_array, amplitude);
xlabel('Cell size (p)')
ylabel('F_0 (k_BT)')