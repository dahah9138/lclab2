clear
clc
clf
%% Parameters
filepath = 'D:/dev/lclab2/data/knot/test';
fileprefix = '';
vortex_type = 'ply';
videoname0 = [filepath '\video'];
rotation_duration = 5;
single_knot_color = 0;
%% Capture settings
% 1 : Only knot
% 2 : Only SB
% 3 : Only preimages
% 4 : Preimages and knot
capture_list = [4];
radius = 0.03;
cpoints = 32;
fps = 20;
scale = 1;
%% End of parameter list

side_deg = 0;
viewangle_side = [side_deg,0];
viewangle_top = [0,90];
viewangle_111 = [45,45];

% Count the number of frames using directory size
a=dir([filepath '/*SB*.ply']);
numframes = size(a,1)/2;
% Frames to stop and rotate on
rotation_frames = [1,numframes];

for capture=capture_list

if capture == 1
    videoname = [videoname0, '_only_knot'];
elseif capture == 2
    videoname = [videoname0, '_sb'];
elseif capture == 3
    videoname = [videoname0, '_preimages'];
elseif capture == 4
    videoname = [videoname0, '_knot_and_preimages'];
end

video_side = VideoWriter([videoname,'_side']);
video_side.FrameRate = fps;
video_side.Quality = 100;
open(video_side);

video_top = VideoWriter([videoname,'_top']);
video_top.FrameRate = fps;
video_top.Quality = 100;
open(video_top);

video_111 = VideoWriter([videoname,'_111']);
video_111.FrameRate = fps;
video_111.Quality = 100;
open(video_111);


for k=1:numframes
    %if k == 331
    %    continue
    %end
    clf('reset')
    %surfaceMeshShow(mesh,Title="Imported Mesh")

    f = gcf;
    f.Position = [100 100 840 800];

    % Read in the meshes
    %draw_mesh([filepath,'/',fileprefix,int2str(k-1),'_vortex.ply'], 0,[1,0,0]);
    % Read in the meshes
    if capture == 0 || capture == 1 || capture == 4
        if single_knot_color==1
            if vortex_type == 'bin'
                f_vortex_knot([filepath,'/',fileprefix,int2str(k-1),'_vortex.bin'],radius,cpoints,0,[1,0,0])
            elseif vortex_type == 'ply'
                draw_mesh_multi([filepath,'/',fileprefix,int2str(k-1),'_vortex'], 0,[0,0,1],alpha);
            end
        else
            if vortex_type == 'bin'
                f_vortex_knot([filepath,'/',fileprefix,int2str(k-1),'_vortex.bin'],radius,cpoints,0)
            elseif vortex_type == 'ply'
                draw_mesh_multi([filepath,'/',fileprefix,int2str(k-1),'_vortex'], 0);
            end
        end
    end
    
    if capture == 0 || capture == 2
        hold on
        draw_mesh([filepath,'/',fileprefix,int2str(k-1),'_SB_pos.ply'], 0,[0,0,1]);
        hold on
        draw_mesh([filepath,'/',fileprefix,int2str(k-1),'_SB_neg.ply'], 0,[1,1,0]);
    end
    
    if capture == 0 || capture == 3 || capture == 4
        hold on
        % Auto detect number of preimages drawn
        draw_preimages([filepath,'/',fileprefix,int2str(k-1),'_preim_'], 0);
    end
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    axis off
    
    set(ax,'CameraViewAngleMode','Manual')

    axis equal
    axis off
    camlight
    camproj('p')

    bd = 1+scale;
    xlim([-bd,bd])
    ylim([-bd,bd])
    zlim([-bd-1,bd])

    if ismember(k,rotation_frames)
        view([side_deg,0])
        dead_frame = getframe(f);
        for stall=1:fps
            writeVideo(video_side,dead_frame);
        end
        for rot=linspace(0,360,rotation_duration*fps)
            view([side_deg+rot,0])
            frame = getframe(f);
            writeVideo(video_side, frame);
        end
        for stall=1:fps
        writeVideo(video_side,dead_frame);
        end
    else
        view(viewangle_side)
        frame = getframe(f);
        writeVideo(video_side, frame);
    end

    view(viewangle_top)
    frame = getframe(f);
    writeVideo(video_top, frame);

    if ismember(k,rotation_frames)
        % Add 1 second dead time before and after rotation
        view([45,45])
        dead_frame = getframe(f);
        for stall=1:fps
            writeVideo(video_111, dead_frame);
        end
        for rot=linspace(0,360,rotation_duration*fps)
            view([45+rot,45])
            frame = getframe(f);
            writeVideo(video_111, frame);
        end
        for stall=1:fps
            writeVideo(video_111, dead_frame);
        end
    else
        view(viewangle_111)
        frame = getframe(f);
        writeVideo(video_111, frame);
    end
end

close(video_side);
close(video_top);
close(video_111);
end

function draw_mesh(file, draw_points, color)
    mesh = readSurfaceMesh(file);
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
                'SpecularStrength',0.3,...
                'AmbientStrength',0.55);
        else
            patch('Vertices', mesh.Vertices, 'Faces',mesh.Faces, 'FaceColor', color,...
            'FaceLighting','flat',...
                'EdgeColor','none',...
                'SpecularStrength',0.3,...
                'AmbientStrength',0.55);
        end
    end
end
function draw_mesh_multi(file, draw_points, color,alpha)
    % Assumes that file path, file prefix, and current frame was passed
    % Loads all components for the current vortex knot
    % file = path/<frame>_vortex
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
