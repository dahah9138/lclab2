% Used to modify .lmt data in matlab and then read it back into .lmt format
% Does not currently work (problem with fwrite_foel...)
clear
clc

file_name = '../data/voltage_dimer2.lmt';

[data, ver] = fread_foel(file_name);

% Things to change:
% Voxels [5]
% celldims [7]
% nn [10]
% voltage [11]

celldims = data{7,2};
nn = data{10,2};
vv = data{11,2};


data{5,2} = size(vv);
data{7,2} = celldims;
data{10,2} = nn;
data{11,2} = vv;

% Save data back to file
fwrite_foel(file_name, data, ver);