function f_vortex_cyl(file,radius,cpoints,color)

%% Parameters
col_specified = 1;

if ~exist('color','var')
    col_specified = 0;
end
%% End of parameters

fid = fopen(file, 'r');
Nc = fread(fid, 1, 'int32');

[xs,ys,zs] = sphere;

for k=1:Nc
    % number of edges in the component
    N = fread(fid, 1,'int32');
    edges = fread(fid, 6*N, 'double');
    edges = reshape(edges, [3 2 N]);
    if col_specified ~= 1
        color = hsv2rgb([k/Nc,1,1]);
    end
    for i=1:N
        % Plot the cylinder
        %[x,y,z]=tubeplot([edges(1,:,i); edges(2,:,i); edges(3,:,i)],radius,cpoints);
        %[c,end1,end2] = Cylinder(edges(:,1,i),edges(:,2,i),radius,cpoints,[1,0,0],0,0);
        %surf(x,y,z,'FaceColor',color,'EdgeColor','none')
        %hold on
        % Plot a sphere at each point
        %surf(radius*xs+edges(1,1,i),radius*ys+edges(2,1,i),radius*zs+edges(3,1,i),'FaceColor',color,'EdgeColor','none')
        %hold on
    end
    plot3(cat(1,squeeze(edges(1,1,:)),squeeze(edges(1,2,end)),squeeze(edges(1,1,1))),...
        cat(1,squeeze(edges(2,1,:)),squeeze(edges(2,2,end)),squeeze(edges(2,1,1))),...
        cat(1,squeeze(edges(3,1,:)),squeeze(edges(3,2,end)),squeeze(edges(3,1,1))))
    hold on
end

fclose(fid);
hold off
end
