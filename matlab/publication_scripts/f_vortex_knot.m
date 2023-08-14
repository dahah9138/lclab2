function components=f_vortex_knot(file,radius,cpoints,pts_only,color,draw)


col_specified = 1;

if ~exist('color','var')
    col_specified = 0;
end

if ~exist('pts_only','var')
    pts_only = 0;
end

if ~exist('draw','var')
    draw = 1;
end

fid = fopen(file, 'r');
Nc = fread(fid, 1, 'int32');

[xs,ys,zs] = sphere;
components = [];

for k=1:Nc
    % number of edges in the component
    N = fread(fid, 1,'int32');
    points = fread(fid, 3*N, 'double');
    
    points = reshape(points, [3 N]);

    components(end+1) = N;
    
    if draw==1
    if col_specified ~= 1
        color = hsv2rgb([k/Nc,1,1]);
    end
    for i=1:N
        % Plot the cylinder
        %[x,y,z]=tubeplot([edges(1,:,i); edges(2,:,i); edges(3,:,i)],radius,cpoints);
        %[c,end1,end2] = Cylinder(edges(:,1,i),edges(:,2,i),radius,cpoints,[1,0,0],0,0);
        %hold on
        % Plot a sphere at each point
        surf(.95*radius*xs+points(1,i),.95*radius*ys+points(2,i),.95*radius*zs+points(3,i),'FaceColor',color,'EdgeColor','none')
        hold on
    end
    if pts_only==0
        % Connect the end
        points = cat(2,points,points(:,1));
        points_interp = interparc(1000,points(1,:),points(2,:),points(3,:),'spline');
    
        [x,y,z]=tubeplot(points_interp',radius,cpoints);
        surf(x,y,z,'FaceColor',color,'EdgeColor','none')
        hold on
    end
    end
end

fclose(fid);
hold off
end
