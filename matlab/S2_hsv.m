% % Color-coded S2 based on orientation

clear
clf
format long e
hold all

% reference vector
cthetad = [];
cphid = [];

% Double angle counting
S2_Z2 = 1;

% viewangle = [0 0];
% viewangle = [36 30];
viewangle = [0,90];

conetheta = cthetad*pi/180;
conephi = cphid*pi/180;

k = 10;
map = hsv(2^k);

% plot axes arrows
arrowlength = 2.5;
fontsize = 36;
stemwidth = 0.03;

% lighting
view(viewangle)



for cone = 1:length(conetheta)
    
% [x1,y1,z1] = cylinder([.25 0],128);
[x1,y1,z1] = cylinder([.25 0],128);

xc = x1;
yc = y1;

zc = z1+.8;

% set up rotation matrix1:
rotationMatrix1 = [cos(conetheta(cone))  0  sin(conetheta(cone));...
                   0                  1  0;...
                   -sin(conetheta(cone)) 0  cos(conetheta(cone))];
% set up rotation matrix2:
rotationMatrix2 = [cos(conephi(cone))  -sin(conephi(cone))  0;...
                   sin(conephi(cone))  cos(conephi(cone))   0;...
                   0                0                 1];
% get points at the two rings and rotate them separately about x:

positionNew1 = rotationMatrix2*rotationMatrix1*[xc(1,:); yc(1,:); zc(1,:)];
positionNew2 = rotationMatrix2*rotationMatrix1*[xc(2,:); yc(2,:); zc(2,:)];

xc = [positionNew1(1,:); positionNew2(1,:)];
yc = [positionNew1(2,:); positionNew2(2,:)];
zc = [positionNew1(3,:); positionNew2(3,:)];

arrow_color = zeros(1,3);
arrow_color(1) = (conephi(cone))/(2*pi); 
if arrow_color(1)<0
    arrow_color(1) = (conephi(cone)+pi)/(2*pi); 
end
arrow_color(2) = conetheta(cone)/(pi/2);
arrow_color(3) = 2 - conetheta(cone)/(pi/2);

for i = 1:3
    if arrow_color(i) > 1
        arrow_color(i) = 1;
    end
end
arrow_color = hsv2rgb(arrow_color);

hold on
s = surf(xc,yc,zc,'EdgeColor','none','FaceColor',arrow_color);
fvc = surf2patch(s);
delete(s)
patch(fvc,...
    'FaceAlpha',1,...
    'FaceColor',arrow_color,...
    'FaceLighting','flat',...
    'EdgeColor','none',...
    'SpecularStrength',0.8,...
    'AmbientStrength',0.55);

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
s1=surf(X,Y,Z,hsv2rgb(col_rs)); % color-code accroding to theta
s1.EdgeColor = 'none';
s1.SpecularStrength = 0.1;
s1.AmbientStrength= 0.4;
s1.FaceAlpha = 1;



colormap(map)
axis equal off
grid off
camlight

fig = gcf;
fig.Color = [1 1 1];
