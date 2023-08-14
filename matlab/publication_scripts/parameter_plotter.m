clear
clc
clf

filepath = 'D:/dev/lclab2/data/knot/test';
fileprefix = '';

% Choose which items to compute
show_coulomb = 1;
show_length = 1;
show_FO_energy = 1;
show_dist = 1;

% Save plots generated in filepath
save_data = 0;

% Find the number of energies to calculate
a=dir([filepath '/*SB*.ply']);
nEnergies = size(a,1)/2;

% Next loop through and compute the energy
coul_energy = zeros(nEnergies,1);
min_distances = zeros(nEnergies,1);
max_distances = zeros(nEnergies,1);
dist_frob_norm = zeros(nEnergies,1);
knot_lengths = cell(nEnergies,1);
knot_lengths_cumul = cell(nEnergies,1);
step = 500;
time = linspace(0,step*nEnergies,nEnergies);
for i = 1:nEnergies
    knot = load_knot([filepath, '/', fileprefix,num2str(i-1),'_vortex.txt']);
    if show_coulomb
        coul_energy(i) = compute_coulomb_energy(knot);
    end
    if show_length
        [knot_lengths{i},knot_lengths_cumul{i}] = compute_piecewise_knot_length(knot);
    end
    if show_dist
        [min_distances(i),max_distances(i),dist_frob_norm(i)] = compute_distance(knot);
    end
end

if show_coulomb
    fig_coulomb = figure(1);
    scatter(time,coul_energy, 'k', 'filled')
    xlabel('Time (a.u.)')
    ylabel('E_C (a.u.)')
    if save_data==1
        savefig(fig_coulomb,[filepath '/coulomb_energy_naive.fig'])
        exportgraphics(fig_coulomb,[filepath '/coulomb_energy_naive.png'],"Resolution", 300)
    end
end

if show_length
    % First create the figure for cumulative length
    fig_length = figure(show_coulomb + 1);
 
    % Find the maximum number of components over all frames
    lmax = 0;
    comp_max = 1;
    for i=1:nEnergies
        flength = cell2mat(knot_lengths(i));
        lmax = max([max(flength),lmax]);
        comp_max = max([length(flength),comp_max]);
    end


    for i=1:nEnergies
        flength = cell2mat(knot_lengths_cumul(i));
        time_temp = linspace(time(i),time(i),length(flength));
        scatter(time_temp,flength,'k','filled')
        hold on
    end
    ylim([0,lmax*1.3])
    xlabel('Time (a.u.)')
    ylabel('Cumulative Length (a.u.)')

    if save_data==1
        savefig(fig_length,[filepath '/cumul_length.fig'])
        exportgraphics(fig_length,[filepath '/cumul_length.png'],"Resolution", 300)
    end

    fig_length = figure(show_coulomb + 2);
    for i=1:nEnergies
        flength = cell2mat(knot_lengths(i));
        time_temp = linspace(time(i),time(i),length(flength));
        % Choose color at frame by number of components
        cnum = length(flength);
        col = hsv2rgb([(cnum-1)/comp_max,1,1]);
        scatter(time_temp,flength,30,col,'filled')
        hold on
    end
    ylim([0,lmax*1.3])
    xlabel('Time (a.u.)')
    ylabel('Length (a.u.)')

    if save_data==1
        savefig(fig_length,[filepath '/length.fig'])
        exportgraphics(fig_length,[filepath '/length.png'],"Resolution", 300)
    end


end

if show_FO_energy
    fig_FO_energy = figure(2*show_length + show_coulomb + 1);
    [FO_energy, npoints] = load_energy(filepath);
    time_full = linspace(0,step*npoints,npoints);
    scatter(time_full,FO_energy, 'filled', 'k')
    xlabel('Time (a.u.)')
    ylabel('Energy (Kp)')
    if save_data==1
        savefig(fig_FO_energy,[filepath '/FO_energy_plot.fig'])
        exportgraphics(fig_FO_energy,[filepath '/FO_energy_plot.png'],"Resolution", 300)
    end
end

if show_dist
    fig_dist = figure(2*show_length + show_coulomb + show_FO_energy + 1);
    scatter(time, min_distances, 'b', 'filled')
    hold on
    scatter(time, max_distances, 'r', 'filled')
    hold on
    scatter(time, dist_frob_norm, 'k', 'filled')
    xlabel('Time (a.u.)')
    ylabel('Distance (p)')
    legend('d_{min}', 'd_{max}', 'd_{Frob}')
    if save_data==1
        savefig(fig_dist,[filepath '/distances.fig'])
        exportgraphics(fig_dist,[filepath '/distances.png'],"Resolution", 300)
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%% End of script
%%%%%%%%%%%%%%%%%%%%%%

% Function that loads knot from a txt file given to it
function knot = load_knot(file_txt)
    data = readlines(file_txt);
    points = zeros(length(data), 3);
    % First component begins at index 1
    comp_split = [0];
    for k=1:length(data)
        if data{k} ~= ""
            points(k,:) = cell2mat(textscan(data{k}, '%f %f %f'));
        else
            comp_split(end+1) = k;
        end
    end
    
    knot = {};
    % Separate components
    for c = 1:(length(comp_split)-1)
        c1 = comp_split(c)+1;
        c2 = comp_split(c+1)-1;
        component = points(c1:c2,:);
        if length(component) > 0
            knot{end+1} = component;
        end
    end

end

function E = compute_coulomb_energy(knot)
    % Merge components together to one point list
    nPoints = 0;
    for C=knot
        nPoints = nPoints + length(cell2mat(C));
    end
    
    R = zeros(nPoints,3);
    k=1;
    for C=knot
        comp = cell2mat(C);
        for pt_i=1:length(comp(:,1))
            R(k,:) = comp(pt_i,:);
            k = k+1;
        end
    end
    
    E = 0;
    for i=2:nPoints
        for j=1:(i-1)
            % Compute the Coulomb energy
            E = E + 1/sqrt(sum((R(i,:)-R(j,:)).^2));
        end
    end
end

function L = compute_total_knot_length(knot)
    % Merge components together to one point list
    nPoints = 0;
    for C=knot
        nPoints = nPoints + length(cell2mat(C));
    end
    
    L = nPoints;
end

function [Li,Li_cumul] = compute_piecewise_knot_length(knot)
    % Merge components together to one point list
    Li = zeros(length(knot),1);
    Li_cumul = zeros(length(knot),1);
    cumulant = 0;
    for C=1:length(knot)
        L0 = length(cell2mat(knot(C)));
        cumulant = cumulant + L0;

        Li_cumul(C) = cumulant;
        Li(C) = L0; 
    end

end

function [energy,npoints] = load_energy(filepath)
    a = dir([filepath '/*energy*.bin']);
    fid = fopen([a.folder '/' a.name]);
    energy = fread(fid, 'double');
    fclose(fid);
    npoints = length(energy);
end

function [dist_min,dist_max, dist_matrix_Frob_norm] = compute_distance(knot)
    nComp = length(knot);
    % Only continue of there are more than one component
    if nComp == 1
        dist_min = 0;
    end

    % Compute the center of mass of each component
    COM_list = zeros(nComp,3);
    for i=1:nComp
        comp = cell2mat(knot(i));
        comp_COM = zeros(1,3);
        comp_len = length(comp(:,1));
        for j=1:comp_len
            comp_COM = comp_COM + comp(j,:);
        end
        COM_list(i,:) = comp_COM/comp_len;
    end
    


    % Find the minimal distance between two components
    if nComp > 1
        dist_min = 10000;
    end
    dist_max = 0;
    dist_matrix_Frob_norm = 0;

    for i=2:nComp
        for j=1:(i-1)
            Ri = COM_list(i,:);
            Rj = COM_list(j,:);
            dij = sqrt(sum((Ri-Rj).^2));
            dist_matrix_Frob_norm = dist_matrix_Frob_norm + 2*dij^2;
            if dij < dist_min && nComp > 1
                dist_min = dij;
            end
            if dij > dist_max
                dist_max = dij;
            end
        end
    end

    dist_matrix_Frob_norm = sqrt(dist_matrix_Frob_norm);
end

