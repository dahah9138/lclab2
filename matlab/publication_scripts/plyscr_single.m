clear
clc
clf
filepath = 'D:/dev/lclab2/data/knot/test';
fileprefix = '';
% Allowed types:
% - ply
% - bin
vortex_type = 'ply';
viewangle = [260,20];
k = 1; % Select which frame to load
use_points = 0;
%% Capture settings
% 0 : Everything
% 1 : Only knot
% 2 : Only SB
% 3 : Only preimages
% 4 : Preimages and knot
capture = 4;
single_knot_color = 0;
radius = 0.03;
cpoints = 32;
alpha = 1;
scale = 1.;
% Draw the vortex knot points on the S2 sphere
spoints = 0;
videoname = 'spheretest';
rotation_duration = 10;
fps = 30;
f = gcf;

set(f,'color','w');

% Read in the meshes
if capture == 0 || capture == 1 || capture == 4
    if single_knot_color==1
        if vortex_type == 'bin'
            f_vortex_knot([filepath,'/',fileprefix,int2str(k-1),'_vortex.bin'],radius,cpoints,0,[1,0,0])
        elseif vortex_type == 'ply'
            draw_mesh_multi([filepath,'/',fileprefix,int2str(k-1),'_vortex'], use_points,[0,0,1],alpha);
        end
    else
        if vortex_type == 'bin'
            f_vortex_knot([filepath,'/',fileprefix,int2str(k-1),'_vortex.bin'],radius,cpoints,0)
        elseif vortex_type == 'ply'
            draw_mesh_multi([filepath,'/',fileprefix,int2str(k-1),'_vortex'], use_points);
        end
    end
end

if capture == 0 || capture == 2
    hold on
    draw_mesh([filepath,'/',fileprefix,int2str(k-1),'_SB_pos.ply'], use_points,[0,0,1],alpha);
    hold on
    draw_mesh([filepath,'/',fileprefix,int2str(k-1),'_SB_neg.ply'], use_points,[1,1,0],alpha);
end

if capture == 0 || capture == 3 || capture == 4
    hold on
    % Auto detect number of preimages drawn
    draw_preimages([filepath,'/',fileprefix,int2str(k-1),'_preim_'], use_points);
end

box on
ax = gca;
ax.BoxStyle = 'full';
axis off

axis equal
axis off
view(viewangle)
camlight
camproj('p')

bd = 1+scale;
xlim([-bd,bd])
ylim([-bd,bd])
zlim([-bd,bd])

if spoints==1
    components=f_vortex_knot([filepath,'/',fileprefix,int2str(k-1),'_vortex.bin'],radius,cpoints,0,[1,0,0],0);
    % Try to load refined vortex line
    fid = fopen([filepath,'/',fileprefix,int2str(k-1),'_sphere.bin'], 'r');
    Npoints = fread(fid, 1, 'int32');
    sphere_data = fread(fid,[2,Npoints],'double');
    fclose(fid);
    theta = sphere_data(1,:);
    phi = sphere_data(2,:)+pi;
    color = zeros(Npoints,3);
    color(:,1) = phi/2/pi;
    color(:,2) = theta/(pi/2);
    color(:,3) = 2 - theta/(pi/2);
    color(color > 1) = 1;
    color = hsv2rgb(color);
    f2 = figure(2);

    % Draw points
    [x,y,z]=sph2cart(phi,pi/2-theta,1);
    %scatter3(x,y,z,[],color,'filled','MarkerEdgeColor','k');
    scatter3sph(x,y,z,'color',color,'size',0.03,'transp',1);
    hold on
    % Draw lines
    cumulant=1;
    ctr=0;
    ncomps = length(components);
    for c=components
        %c
        c_upper=cumulant+c-1;
        phi_sub = phi(cumulant:c_upper);
        theta_sub = theta(cumulant:c_upper);
        [xc,yc,zc]=sph2cart(phi_sub,pi/2-theta_sub,1);
        xc(end+1) = xc(1); yc(end+1) = yc(1); zc(end+1) = zc(1);
        plot3(xc,yc,zc,'color',hsv2rgb([ctr/ncomps 1 1]));
        hold on
        cumulant = c_upper+1;
        ctr = ctr + 1;
    end

    draw_S2_sphere(.98)


    box on
    ax = gca;
    ax.BoxStyle = 'full';
    axis off
    set(ax,'CameraViewAngleMode','Manual')
    
    set(f2,'color','w');
    f2.Position = [100 100 840 800];
    axis equal
    axis off
    view(viewangle)
    camlight
    camproj('p')
    xlim([-1,1])
    ylim([-1,1])
    zlim([-1,1])

    video = VideoWriter([videoname,'_side'],'MPEG-4');
    video.FrameRate = fps;
    open(video);

    for rot=linspace(0,360,rotation_duration*fps)
        view([rot,45])
        frame = getframe(f2);
        writeVideo(video, frame);
    end
    close(video);

end

function draw_mesh(file, draw_points, color,alpha)
    mesh = readSurfaceMesh(file);

    if ~exist('alpha','var')
        alpha = 1;
    end

    if draw_points==1
        if exist('color','var')
            s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled','MarkerFaceColor',color);
        else
            s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled');
        end
        s.SizeData = 4;
    else
        if ~exist('color','var')
            patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceVertexCData', mesh.VertexColors,...
            'FaceColor','interp',...
            'FaceLighting','flat',...
                'EdgeColor','none',...
                'FaceAlpha',alpha,...
                'SpecularStrength',0.3,...
                'AmbientStrength',0.55);
        else
            patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceColor', color,...
            'FaceLighting','flat',...
                'EdgeColor','none',...
                'FaceAlpha',alpha,...
                'SpecularStrength',0.3,...
                'AmbientStrength',0.55);

        hold on
        end
    end
end

function [S2points,scols] = draw_mesh_multi(file, draw_points, color,alpha)
    % Assumes that file path, file prefix, and current frame was passed
    % Loads all components for the current vortex knot
    % file = path/<frame>_vortex
    a=dir([file '*_of_*.ply']);
    num_components = size(a,1);
    S2points = 0;
    scols = 0;

    for c=1:num_components
        fname = [file num2str(c) '_of_', num2str(num_components) '.ply'];
        mesh = readSurfaceMesh(fname);

        if ~exist('alpha','var')
            alpha = 1;
        end
    
        if draw_points==1
            if exist('color','var')
                s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled','MarkerFaceColor',color);
            else
                s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled');
            end
            s.SizeData = 4;
        else
            if ~exist('color','var')
                patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceVertexCData', mesh.VertexColors,...
                'FaceColor','interp',...
                'FaceLighting','flat',...
                    'EdgeColor','none',...
                    'FaceAlpha',alpha,...
                    'SpecularStrength',0.3,...
                    'AmbientStrength',0.55);
                sc = mesh.VertexColors;
                sp = rgb2hsv(sc);
                
                % Get theta comp
                npts = length(sp(:,1));
                theta = zeros(1,npts);
                for i=1:npts
                    if sp(i,2) <= sp(i,3)
                        theta(i)=(pi/2) * sp(i,2);
                    else
                        theta(i)=pi-pi/2*sp(i,3);
                    end
                end
                %theta1=(pi/2) * sp(:,2);
                %theta2=pi - pi/2*sp(:,3); % attach component
                %theta= min(cat(1,theta1',theta2'))';
                % Convert hsv to theta phi
                sp(:,1) = 2*pi*sp(:,1);
                sp(:,2) = theta;
                sp(:,3) = c; % attach component
                if S2points==0
                    S2points = sp;
                    scols=sc;
                else
                    S2points = cat(1,S2points, sp);
                    scols = cat(1,scols, sc);
                end
                
            else
                % Make a color wheel color
                theta = (c-1)/num_components;
                color = hsv2rgb([theta,1,1]);
                patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceColor', color,...
                'FaceLighting','flat',...
                    'EdgeColor','none',...
                    'FaceAlpha',alpha,...
                    'SpecularStrength',0.3,...
                    'AmbientStrength',0.55);
            end
        end
    end
end

function draw_preimages(file, draw_points, color,alpha)
    % Assumes that file path, file prefix, and current frame was passed
    % Loads all components for the current vortex knot
    % file = path/<prefix><frame>_preim
    a=dir([file '*_of_*.ply']);
    num_components = size(a,1);

    for c=1:num_components
        fname = [file num2str(c) '_of_', num2str(num_components) '.ply'];
        mesh = readSurfaceMesh(fname);

        if ~exist('alpha','var')
            alpha = 1;
        end
    
        if draw_points==1
            if exist('color','var')
                s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled','MarkerFaceColor',color);
            else
                s=scatter3(mesh.Vertices(:,1),mesh.Vertices(:,2),mesh.Vertices(:,3),'filled');
            end
            s.SizeData = 4;
        else
            if ~exist('color','var')
                
                patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceVertexCData', mesh.VertexColors,...
                'FaceColor','interp',...
                'FaceLighting','flat',...
                    'EdgeColor','none',...
                    'FaceAlpha',alpha,...
                    'SpecularStrength',0.3,...
                    'AmbientStrength',0.55);
                
            else
                % Make a color wheel color
                theta = (c-1)/num_components;
                color = hsv2rgb([theta,1,1]);
                patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceColor', color,...
                'FaceLighting','flat',...
                    'EdgeColor','none',...
                    'FaceAlpha',alpha,...
                    'SpecularStrength',0.3,...
                    'AmbientStrength',0.55);
            end
        end
    end
end

function draw_S2_sphere(r,S2_Z2,alpha)
    if ~exist('S2_Z2','var')
        S2_Z2 = 0;
    end
    if ~exist('alpha','var')
        alpha = 0.5;
    end

    [X,Y,Z]=sphere(500);
    [m,n] = size(X);
    [azi,ele,~] = cart2sph(X,Y,Z);
    % Z(Z >= cos(bound_d*pi/180)) = nan; % remove the region with theta<bound
    
    azi_rs = reshape(azi,[m*n,1]);
    azi_rs(azi_rs<0) = azi_rs(azi_rs<0)+2*pi;
    the_rs = pi-(reshape(ele,[m*n,1])+pi/2);
    col = zeros(m*n,3);
    
    % hsv: s = 0 at the north pole; v = 0 at the south pole;
    for i = 1:m*n
        
        if S2_Z2
            col(i,1) = mod(2*azi_rs(i),2*pi)/(2*pi);
            col(i,2) = 1;
            col(i,3) = 1;
        else
            col(i,1) = (azi_rs(i))/(2*pi);
            col(i,2) = the_rs(i)/(pi/2);
            col(i,3) = 2 - the_rs(i)/(pi/2);
        end
    
    end
    col(col > 1) = 1;
    col_rs = reshape(col,[m,n,3]);
    % figure
    s1=surf(X*r,Y*r,Z*r,hsv2rgb(col_rs)); % color-code accroding to theta
    s1.EdgeColor = 'none';
    s1.SpecularStrength = 0.1;
    s1.AmbientStrength= 0.4;
    s1.FaceAlpha = alpha;
    k = 10;
    map = hsv(2^k);
    colormap(map)
end