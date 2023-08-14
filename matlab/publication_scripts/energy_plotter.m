clear
clc
clf
files = {};

%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%
save_plot = 0;
step = 500;

% Main file path where plot will be saved
primary_filepath = 'D:/dev/lclab2/data/knot/test';

% Add as many files as desired here
files{end+1} = primary_filepath;
%%%%%%%%%%%%%%

f=figure(1);
for file=files
    [energy, npoints] = load_energy(file{1},step);
    
    X = linspace(0,step*npoints,npoints);
    scatter(X,energy, 'filled', 'k')
    hold on
end
xlabel('Time (a.u.)')
ylabel('Energy (Kp)')
%legend('E = 0.6','E = 0.7','E = 0.8')
hold off

if save_plot==1
    savefig(f,[primary_filepath '/energy_plot.fig'])
    exportgraphics(f,[primary_filepath '/energy_plot.png'],"Resolution", 300)
end

function [energy,npoints] = load_energy(filepath,step)
    a = dir([filepath '/*energy*.bin']);
    fid = fopen([a.folder '/' a.name]);
    energy = fread(fid, 'double');
    fclose(fid);
    npoints = length(energy);
end