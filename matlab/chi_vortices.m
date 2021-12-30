% Plotter for preimage, xsections
clear
clc

tic
% figure
clf
hold all
delete(findall(gcf,'Type','light'))

file = 'dimer';
file_chi = strcat(file, '_chi');

Starget = 0.533;
epsilon = 0.1;

% Looks the same either way for heliknoton
three_fields = 0;

% Sampling rate for the xz-xsection
ratexs = 4; 

% Sampling rate for isosurfaces of preimages
sample = 5; 
xs_xz = 0;
xs_xy = 0;
xs_yz = 0;

load(file, 'nn');
load(file_chi, 'chi');

% make nn same size as chi
nn = nn(2:end-1,2:end-1,2:end-1, :);

if three_fields

    Starget = 0.0;
    epsilon = 0.2;
    
    % compute tau

    tau = cross(nn, chi, 4);
    modtau = sqrt(tau(:,:,:,1).^2+tau(:,:,:,2).^2+tau(:,:,:,3).^2);

    tau(:,:,:,1) = tau(:,:,:,1)./modtau;
    tau(:,:,:,2) = tau(:,:,:,2)./modtau;
    tau(:,:,:,3) = tau(:,:,:,3)./modtau;

    [m0,n0,p0,~] = size(nn);

    Sn = scalarOP(nn);
    Sc = scalarOP(chi);
    St = scalarOP(tau);
    
    S1 = 0.5*(Sc + St) - Sn;
else
    Sc = scalarOP(chi);
    S1 = Sc;
end

clear Sc St Sn chi tau

viewangle = [0,90];
xsection_alpha = 0.4;

preimg_alpha = 1;


[m,n,p,d] = size(nn);
m_int = m*sample;
n_int = n*sample;
p_int = p*sample;

xx = linspace(1,m,m*sample);
yy = linspace(1,n,n*sample);
zz = linspace(1,p,p*sample);

[xx,yy,zz] = ndgrid(xx,yy,zz);

disp('Created ndgrid')

nn_int(:,:,:,1) = interp3(nn(:,:,:,1),yy,xx,zz,'cubic');
nn_int(:,:,:,2) = interp3(nn(:,:,:,2),yy,xx,zz,'cubic');
nn_int(:,:,:,3) = interp3(nn(:,:,:,3),yy,xx,zz,'cubic');


s1_int(:,:,:) = interp3(S1(:,:,:),yy,xx,zz,'cubic');


disp('interpolated')
% % normalize director field after interpolation
modnn = sqrt(nn_int(:,:,:,1).^2+nn_int(:,:,:,2).^2+nn_int(:,:,:,3).^2);
nn_int(:,:,:,1) = nn_int(:,:,:,1)./modnn;
nn_int(:,:,:,2) = nn_int(:,:,:,2)./modnn;
nn_int(:,:,:,3) = nn_int(:,:,:,3)./modnn;

disp('finished normalization')

nn2 = nn_int;
[ma,na,pa,~] = size(nn2);


diffmag = abs(s1_int-Starget);

disp('Finished diffmag')

fv = isosurface(xx,yy,zz,diffmag,epsilon);

disp('Created isosurface')


if numel(fv.vertices)>3

% % old color scheme
%     color = map(end-round(abs(theta(cone))/pi*(2^k-1)),:);

% % HSV color scheme
% %     linear
    color = zeros(1,3);
    color(1) = 0;
    color(2) = 1/2;
    color(3) = 1;
    color(color > 1) = 1;
    color = hsv2rgb(color);

patch(fv,...
        'FaceColor',color,...
        'FaceLighting','gouraud',...
        'FaceAlpha',preimg_alpha,...
        'EdgeColor','none',...
        'SpecularStrength',0.4,...
        'AmbientStrength',0.75);



end


% Call f_preimage now
f_preimage(nn, sample);
    


switch xs_xz
    case 1

[yy0,xx0,zz0] = meshgrid(1:n,1:m,1:p);
[yyq,xxq,zzq] = meshgrid((1+n)/2,1:m,1:p);

% % field on the midplane
nnxzmid(:,:,:,1) = interp3(yy0,xx0,zz0,nn(:,:,:,1),yyq,xxq,zzq,'cubic');
nnxzmid(:,:,:,2) = interp3(yy0,xx0,zz0,nn(:,:,:,2),yyq,xxq,zzq,'cubic');
nnxzmid(:,:,:,3) = interp3(yy0,xx0,zz0,nn(:,:,:,3),yyq,xxq,zzq,'cubic');

% % 
x = linspace(1,m,m*ratexs);
z = linspace(1,p,p*ratexs);
% 
[zz,xx] = meshgrid(z,x);
nn_int2(:,:,1) = interp2(squeeze(nnxzmid(:,:,:,1)),zz,xx,'cubic');
nn_int2(:,:,2) = interp2(squeeze(nnxzmid(:,:,:,2)),zz,xx,'cubic');
nn_int2(:,:,3) = interp2(squeeze(nnxzmid(:,:,:,3)),zz,xx,'cubic');

modnn = sqrt(nn_int2(:,:,1).^2+nn_int2(:,:,2).^2+nn_int2(:,:,3).^2);
nn_int2(:,:,1) = nn_int2(:,:,1)./modnn;
nn_int2(:,:,2) = nn_int2(:,:,2)./modnn;
nn_int2(:,:,3) = nn_int2(:,:,3)./modnn;

thetaxs = real(acos(nn_int2(:,:,3))); % theta 
phixs = fatan(nn_int2(:,:,1),nn_int2(:,:,2)); % branch cut at pi, -pi
phixs(isnan(phixs)) = 0;
phixs(phixs < 0) = phixs(phixs < 0) + 2*pi;
CData = zeros(m*ratexs,p*ratexs,3);

CData(:,:,1) = phixs/2/pi;
CData(:,:,2) = thetaxs/(pi/2);
CData(:,:,3) = 2 - thetaxs/(pi/2);
CData(CData > 1) = 1;
CData = hsv2rgb(CData);

S = surf(xx,(n+1)/2*ones(size(xx)),zz,CData);
S.EdgeColor = 'none';
S.AmbientStrength = 0.3;
S.FaceAlpha = xsection_alpha;
    case 0
end

switch xs_xy
    case 1

[yy0,xx0,zz0] = meshgrid(1:n,1:m,1:p);
[yyq,xxq,zzq] = meshgrid(1:n,1:m,(p+1)/2);

% % field on the midplane
nnxymid(:,:,:,1) = interp3(yy0,xx0,zz0,nn(:,:,:,1),yyq,xxq,zzq,'cubic');
nnxymid(:,:,:,2) = interp3(yy0,xx0,zz0,nn(:,:,:,2),yyq,xxq,zzq,'cubic');
nnxymid(:,:,:,3) = interp3(yy0,xx0,zz0,nn(:,:,:,3),yyq,xxq,zzq,'cubic');

% % 
x = linspace(1,m,m*ratexs);
y = linspace(1,n,n*ratexs);
% 
[yy,xx] = meshgrid(y,x);
clear nn_int2
nn_int2(:,:,1) = interp2(squeeze(nnxymid(:,:,:,1)),yy,xx,'cubic');
nn_int2(:,:,2) = interp2(squeeze(nnxymid(:,:,:,2)),yy,xx,'cubic');
nn_int2(:,:,3) = interp2(squeeze(nnxymid(:,:,:,3)),yy,xx,'cubic');

modnn = sqrt(nn_int2(:,:,1).^2+nn_int2(:,:,2).^2+nn_int2(:,:,3).^2);
nn_int2(:,:,1) = nn_int2(:,:,1)./modnn;
nn_int2(:,:,2) = nn_int2(:,:,2)./modnn;
nn_int2(:,:,3) = nn_int2(:,:,3)./modnn;

thetaxs = real(acos(nn_int2(:,:,3))); % theta 
phixs = fatan(nn_int2(:,:,1),nn_int2(:,:,2)); % branch cut at pi, -pi
phixs(isnan(phixs)) = 0;
phixs(phixs < 0) = phixs(phixs < 0) + 2*pi;
CData = zeros(m*ratexs,n*ratexs,3);

CData(:,:,1) = phixs/2/pi;
CData(:,:,2) = thetaxs/(pi/2);
CData(:,:,3) = 2 - thetaxs/(pi/2);
CData(CData > 1) = 1;
CData = hsv2rgb(CData);

S = surf(xx,yy,(p+1)/2*ones(size(xx)),CData);
S.EdgeColor = 'none';
S.AmbientStrength = 0.3;
S.FaceAlpha = xsection_alpha;
    case 0
end

switch xs_yz
    case 1

[yy0,xx0,zz0] = meshgrid(1:n,1:m,1:p);
[yyq,xxq,zzq] = meshgrid(1:n,(m+1)/2,1:p);

% % field on the midplane
nnyzmid(:,:,:,1) = interp3(yy0,xx0,zz0,nn(:,:,:,1),yyq,xxq,zzq,'cubic');
nnyzmid(:,:,:,2) = interp3(yy0,xx0,zz0,nn(:,:,:,2),yyq,xxq,zzq,'cubic');
nnyzmid(:,:,:,3) = interp3(yy0,xx0,zz0,nn(:,:,:,3),yyq,xxq,zzq,'cubic');

% % 
y = linspace(1,n,n*ratexs);
z = linspace(1,p,p*ratexs);
% 
[zz,yy] = meshgrid(z,y);
clear nn_int2
nn_int2(:,:,1) = interp2(squeeze(nnyzmid(:,:,:,1)),zz,yy,'cubic');
nn_int2(:,:,2) = interp2(squeeze(nnyzmid(:,:,:,2)),zz,yy,'cubic');
nn_int2(:,:,3) = interp2(squeeze(nnyzmid(:,:,:,3)),zz,yy,'cubic');

modnn = sqrt(nn_int2(:,:,1).^2+nn_int2(:,:,2).^2+nn_int2(:,:,3).^2);
nn_int2(:,:,1) = nn_int2(:,:,1)./modnn;
nn_int2(:,:,2) = nn_int2(:,:,2)./modnn;
nn_int2(:,:,3) = nn_int2(:,:,3)./modnn;

thetaxs = real(acos(nn_int2(:,:,3))); % theta 
phixs = fatan(nn_int2(:,:,1),nn_int2(:,:,2)); % branch cut at pi, -pi
phixs(isnan(phixs)) = 0;
phixs(phixs < 0) = phixs(phixs < 0) + 2*pi;
CData = zeros(n*ratexs,p*ratexs,3);
CData(:,:,1) = phixs/2/pi;
CData(:,:,2) = thetaxs/(pi/2);
CData(:,:,3) = 2 - thetaxs/(pi/2);
CData(CData > 1) = 1;
CData = hsv2rgb(CData);

S = surf((m+1)/2*ones(size(yy)),yy,zz,CData);
S.EdgeColor = 'none';
S.AmbientStrength = 0.3;
S.FaceAlpha = xsection_alpha;
    case 0
end
 
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

% Assumes sample uses PBCs
function S = scalarOP(n)

nb = 6; % number of neighbors

avgcos2 = (dot(n,n([2:end 1],:,:,:),4).^2+...
           dot(n,n([end 1:end-1],:,:,:),4).^2+...
           dot(n,n(:,[end 1:end-1],:,:),4).^2+...
           dot(n,n(:,[2:end 1],:,:),4).^2+...
           dot(n,n(:,:,[end 1:end-1],:),4).^2+...
           dot(n,n(:,:,[2:end 1],:),4).^2)/nb;
       
S = 0.5*(3*avgcos2-1);
end