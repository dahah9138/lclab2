% Plotter for preimage, xsections
clear 
tic
% figure
clf
hold all
delete(findall(gcf,'Type','light'))

file = 'Q2_sim';

% Cone error
eps = 0.2;

% Sampling rate for isosurfaces of preimages
sample = 2; 
preimg_alpha = 1;

file = strcat(file, '.mat');

try
    load(file,'vv');
    disp('loaded vv')
catch    
    disp('failed to load vv')
end

% Cut off some of nn

[m0,n0,p0,~] = size(vv);

p01 = int32(p0 *3/9);
p02 = int32(p0 *6/9);

cutvv = vv(:,:,p01:p02,:);
%nn = cutnn;


xs_xz = 0;
xs_xy = 0;
xs_yz = 0;
viewangle = [0,90];
xsection_alpha = 0.4;

% Sampling rate for the xz-xsection
ratexs = 2;

[m,n,p,d] = size(vv);
m_int = m*sample;
n_int = n*sample;
p_int = p*sample;

xx = linspace(1,m,m*sample);
yy = linspace(1,n,n*sample);
zz = linspace(1,p,p*sample);

[X,Y,Z] = ndgrid(xx,yy,zz);

disp('Created ndgrid')

vv_int(:,:,:) = interp3(vv(:,:,:),Y,X,Z,'cubic');

disp('interpolated')

phalf = cast(p_int/2,'int32');
vplane = vv_int(:,:,phalf);

[yy,xx] = meshgrid(yy,xx);

S = surf(xx,yy,(p+1)/2*ones(size(xx)), vplane, 'edgecolor','none');
 
% xlabel('x')
% ylabel('y')
% zlabel('z')


box
ax = gca;
ax.BoxStyle = 'full';
axis off

set(gcf,'Color','w')
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
axis equal
xlim([1 m])
ylim([1 n])
zlim([1 p])
view(viewangle)
camlight
camproj('p')
% %%
% %% axis arrows
hold all
axis off
box on
xlim([-2 m+1])
ylim([-2 n+1])
zlim([-2 p+1])
linewidth = 1;


arrow_length = 16;
stemwidth = 0.5;
mArrow3([1 1 1],[arrow_length 1 1],'facealpha', 1, 'stemWidth', stemwidth);
% text(24,2,2,'x','FontSize',12)
mArrow3([1 1 1],[1 arrow_length 1],'facealpha', 1, 'stemWidth', stemwidth);
% text(2,24,2,'y','FontSize',12)
mArrow3([1 1 1],[1 1 arrow_length],'facealpha', 1, 'stemWidth', stemwidth);
% text(2,2,24,'z','FontSize',12)
set(gcf,'Color',1*[1 1 1])

%% axis arrows and box

hold all
scale = 1.0;
small_bd = (m+1)/2-(m+1)/2*scale;
large_bd = (m+1)/2+(m+1)/2*scale;
xlim([small_bd-2 large_bd+1])
ylim([small_bd-2 large_bd+1])
zlim([-2 p+1])
linewidth = 1;


function A=fatan(X,Y)
% fatan.m
% proper arctangent!
% daa 081020
a=(X<0);  % 0 in I, IV; 1 in II, III
b=(-1).^a;  % 1 in I, IV; -1 in II, III
c=(-1).^(Y<0);  % 1 in I, II; -1 in III, IV
A=-pi*a.*b.*c+atan(Y./X);
end
