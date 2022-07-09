% % Finds the path along the max/min of a 3d scalar field
% % e.g., finds the path along certain preimage

clear
format short
clc
clf
hold all

% % load file

file = 'mag_heliknoton_191125f_3pitch_216_chi';
load('mag_heliknoton_191125f_3pitch_216_chi','data')
chidiff = data.chidiff;
% % 

% % smooth data
filter = 'box';
% filter = 'gaussian';
boxsize = 1*[1 1 1];
chidiff = smooth3(chidiff,filter,boxsize);
% % 

% % interpolate
sample = 1;
[m,n,p] = size(chidiff);
[xx,yy,zz] = ndgrid(linspace(1,m,m*sample),...
    linspace(1,n,n*sample),linspace(1,p,p*sample));
chidiff = interp3(chidiff,yy,xx,zz,'cubic');
% % 

% % Make interpolant of the scalar field
F = griddedInterpolant(xx,yy,zz,chidiff);
% % 

% % plot isosurface/preimage
fv = isosurface(xx,yy,zz,chidiff,60);
f = patch(fv,...
        'FaceColor',[1 1 1],...
        'FaceLighting','gouraud',...
        'FaceAlpha',1,...
        'EdgeColor','none',...
        'SpecularStrength',0.5,...
        'AmbientStrength',0.2);
% % 
    
view(0,0)
camlight
axis equal

%% curve.three
% % Find the first point
clf
hold all

% % Give an initial tangent and inital starting point
tangent = [2 1 0];
x0 = [108 108 151.3];
% tangent = tangent./sqrt(tangent(:,1).^2+tangent(:,2).^2+tangent(:,3).^2)
% % 

% % make basis vectors orthonormal to tangent
p(1) = abs(tangent(2));
p(2) = abs(tangent(1));
p(3) = sum(p,2)==0;
perp1 = cross(p,tangent);
perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
perp2 = cross(perp1,tangent);
% % 

% % some parameters
step = 2;
n = 30; % points in each circle
tc = linspace(0,2*pi,n);
r = linspace(0,5,11); % radius
% % 

x = [];
y = [];
z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

    [~,index] = max(F(x,y,z));    
    x0 = [x(index) y(index) z(index)];
    
% %%
xyz = [x0(1); x0(2); x0(3)];

% % Iteratively find subsequent pts
for i = 1:242
    
    step = 0.8;

    p(1) = abs(tangent(2));
    p(2) = abs(tangent(1));
    p(3) = sum(p,2)==0;

    perp1 = cross(p,tangent);
    perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
    perp2 = cross(perp1,tangent);

    n = 50; % points in each circle
    tc = linspace(0,2*pi,n);
    r = linspace(0,3,20); % radius

    x = [];
    y = [];
    z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

%     plot3(x,y,z)

    [~,index] = max(F(x,y,z));

    plot3(x(index),y(index),z(index),'d')
    
    x00 = x0;
    x0 = [x(index) y(index) z(index)];
    tangent = x0 - x00;
    tangent = tangent./sqrt(tangent(1).^2+tangent(2).^2+tangent(3).^2);
    
    xyz = [xyz x0'];

end

% % Smooth curve
xyz = xyz';
xyz = smoothdata(xyz,'gaussian',20);

% % interpolate to finer pts
xarray = interp1(1:size(xyz,1),xyz(:,1),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
yarray = interp1(1:size(xyz,1),xyz(:,2),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
zarray = interp1(1:size(xyz,1),xyz(:,3),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
xyz = [xarray' yarray' zarray'];

% % create tube around the curve
parametric_tube(xyz,1,0.5*[1 1 1],1)
axis equal tight
view(180,0)
camlight
%% curve.two
% % Find the first point
clf
hold all

% % Give an initial tangent and inital starting point
tangent = [0 1 -1];
x0 = [52.43 108 108];
% tangent = tangent./sqrt(tangent(:,1).^2+tangent(:,2).^2+tangent(:,3).^2)
% % 

% % make basis vectors orthonormal to tangent
p(1) = abs(tangent(2));
p(2) = abs(tangent(1));
p(3) = sum(p,2)==0;
perp1 = cross(p,tangent);
perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
perp2 = cross(perp1,tangent);
% % 

% % some parameters
step = 2;
n = 30; % points in each circle
tc = linspace(0,2*pi,n);
r = linspace(0,5,11); % radius
% % 

x = [];
y = [];
z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

    [~,index] = max(F(x,y,z));    
    x0 = [x(index) y(index) z(index)];
    
% %%
xyz = [x0(1); x0(2); x0(3)];

% % Iteratively find subsequent pts
for i = 1:239
    
    step = 0.8;

    p(1) = abs(tangent(2));
    p(2) = abs(tangent(1));
    p(3) = sum(p,2)==0;

    perp1 = cross(p,tangent);
    perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
    perp2 = cross(perp1,tangent);

    n = 50; % points in each circle
    tc = linspace(0,2*pi,n);
    r = linspace(0,3,20); % radius

    x = [];
    y = [];
    z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

%     plot3(x,y,z)

    [~,index] = max(F(x,y,z));

    plot3(x(index),y(index),z(index),'d')
    
    x00 = x0;
    x0 = [x(index) y(index) z(index)];
    tangent = x0 - x00;
    tangent = tangent./sqrt(tangent(1).^2+tangent(2).^2+tangent(3).^2);
    
    xyz = [xyz x0'];

end

% % Smooth curve
xyz = xyz';
xyz = smoothdata(xyz,'gaussian',35);

% % interpolate to finer pts
xarray = interp1(1:size(xyz,1),xyz(:,1),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
yarray = interp1(1:size(xyz,1),xyz(:,2),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
zarray = interp1(1:size(xyz,1),xyz(:,3),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
xyz = [xarray' yarray' zarray'];

% % create tube around the curve
parametric_tube(xyz,1,0.5*[1 1 1],1)
axis equal
view(90,90)
camlight
%% curve.one
% % Find the first point
clf
hold all

% % Give an initial tangent and inital starting point
tangent = [0 1 1];
x0 = [163.4 109 108];
% tangent = tangent./sqrt(tangent(:,1).^2+tangent(:,2).^2+tangent(:,3).^2)
% % 

% % make basis vectors orthonormal to tangent
p(1) = abs(tangent(2));
p(2) = abs(tangent(1));
p(3) = sum(p,2)==0;
perp1 = cross(p,tangent);
perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
perp2 = cross(perp1,tangent);
% % 

% % some parameters
step = 2;
n = 30; % points in each circle
tc = linspace(0,2*pi,n);
r = linspace(0,5,11); % radius
% % 

x = [];
y = [];
z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

    [~,index] = max(F(x,y,z));    
    x0 = [x(index) y(index) z(index)];
    
% %%
xyz = [x0(1); x0(2); x0(3)];

% % Iteratively find subsequent pts
for i = 1:236
    
    step = 0.8;

    p(1) = abs(tangent(2));
    p(2) = abs(tangent(1));
    p(3) = sum(p,2)==0;

    perp1 = cross(p,tangent);
    perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
    perp2 = cross(perp1,tangent);

    n = 50; % points in each circle
    tc = linspace(0,2*pi,n);
    r = linspace(0,3,20); % radius

    x = [];
    y = [];
    z = [];
    for i = 1:length(r)
        x = [x x0(1) + step*tangent(1) + r(i)*(cos(tc)*perp1(1) + sin(tc)*perp2(1))];
        y = [y x0(2) + step*tangent(2) + r(i)*(cos(tc)*perp1(2) + sin(tc)*perp2(2))];
        z = [z x0(3) + step*tangent(3) + r(i)*(cos(tc)*perp1(3) + sin(tc)*perp2(3))];
    end

%     plot3(x,y,z)

    [~,index] = max(F(x,y,z));

%     plot3(x(index),y(index),z(index),'d')
    
    x00 = x0;
    x0 = [x(index) y(index) z(index)];
    tangent = x0 - x00;
    tangent = tangent./sqrt(tangent(1).^2+tangent(2).^2+tangent(3).^2);
    
    xyz = [xyz x0'];

end

% % Smooth curve
xyz = xyz';
xyz = smoothdata(xyz,'gaussian',35);

% % interpolate to finer pts
xarray = interp1(1:size(xyz,1),xyz(:,1),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
yarray = interp1(1:size(xyz,1),xyz(:,2),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
zarray = interp1(1:size(xyz,1),xyz(:,3),linspace(1,size(xyz,1),4*size(xyz,1)),'pchip');
xyz = [xarray' yarray' zarray'];

% % create tube around the curve
parametric_tube(xyz,1,0.5*[1 1 1],1)
axis equal
view(90,90)
camlight

%%

function parametric_tube(xyz,r,color,facealpha)
 
n = 50; % points in each circle
tc = linspace(0,2*pi,n)';
 
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
 
tangent = cat(1,[diff(x),diff(y),diff(z)],nan(1,3));
tangent = tangent./sqrt(tangent(:,1).^2+tangent(:,2).^2+tangent(:,3).^2);
 
% make a normalized vector basis
p(:,3) = abs(tangent(:,1))<1;
p(:,1) = abs(tangent(:,2))>1;
p(:,2) = sum(p,2)==0;
perp1 = cross(p,tangent);
perp1 = perp1./sqrt(perp1(:,1).^2+perp1(:,2).^2+perp1(:,3).^2);
perp2 = cross(perp1,tangent);
 
c(:,:,1) = x+r*cos(tc)'.*perp1(:,1)+r*sin(tc)'.*perp2(:,1);
c(:,:,2) = y+r*cos(tc)'.*perp1(:,2)+r*sin(tc)'.*perp2(:,2);
c(:,:,3) = z+r*cos(tc)'.*perp1(:,3)+r*sin(tc)'.*perp2(:,3);
c(end,:,:) = c(1,:,:);
 

surf(c(:,:,1),c(:,:,2),c(:,:,3),'LineStyle','none','FaceColor',color,'SpecularStrength',1,'FaceAlpha',facealpha,'AmbientStrength',0.3)


end