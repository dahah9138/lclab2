clear
clc
clf

%% Parameters
file = '../data/knot/test/points.bin';
radius = 0.03;
cpoints = 32;
routine = 1;
color = [1,0,0];
%% End of parameters

fid = fopen(file, 'r');
Nc = fread(fid, 1, 'int32');

[xs,ys,zs] = sphere;

for k=1:Nc
    % number of edges in the component
    N = fread(fid, 1,'int32');
    edges = fread(fid, 6*N, 'double');
    edges = reshape(edges, [3 2 N]);
    color = hsv2rgb([k/Nc,1,1]);
    for i=1:N
        % Plot the cylinder
        nh = edges(:,2,i)-edges(:,1,i);
        midpoint = 0.5*(edges(:,2,i)+edges(:,1,i));
        nh_len = dot(nh,nh,1);
        nh = nh/nh_len;
        [x,y,z]=tubeplot([edges(1,:,i); edges(2,:,i); edges(3,:,i)],radius,cpoints);
        %[c,end1,end2] = Cylinder(edges(:,1,i),edges(:,2,i),radius,cpoints,[1,0,0],0,0);
        surf(x,y,z,'FaceColor',color,'EdgeColor','none')
        hold on
    
        % Plot a sphere at each point
        surf(radius*xs+edges(1,1,i),radius*ys+edges(2,1,i),radius*zs+edges(3,1,i),'FaceColor',color,'EdgeColor','none')
        hold on
    end
end

fclose(fid);

hold off
daspect([1,1,1]); camlight;
xlim([-3,3])
ylim([-3,3])
zlim([-3,3])