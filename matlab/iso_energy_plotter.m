% Plotter for preimage, xsections
clear 
tic
% figure
clf
hold all
delete(findall(gcf,'Type','light'))
preimg_alpha = 0.8;
xsection_alpha = 0.6;
viewangle = [0,90];

xs_xy = 0;
xs_xz = 0;
xs_yz = 0;
ratexs = 3;
sample = 4;
use_zeropoint = 0;
use_k22_only = 0;
energy_iso = 1;

% relative file location to iso_energy_plotter.m
file_loc = "";
% file name
file = "Q1-one-const";
% file extension
file_type = ".mat";

% % 5CB % % % ---------------------------------------------------
%k11 = 6.4/6.47;    % elastic constants scaled by K22
%k22 = 3/6.47;
%k33 = 10/6.47;
%k24 = 0;

% One-const free en
k11 = 1;
k22 = 1;
k33 = 1;
k24 = 1;

chirality = 1.0;


try
    load(strcat(file_loc, file, file_type), 'cdims');
catch
    warning("Failed to load cdims")
end

if ~exist('cdims','var')
    
    cdims = zeros(3);
    
    cdims(1) = input("Input x cell dimension: ");
    cdims(2) = input("Input y cell dimension: ");
    cdims(3) = input("Input z cell dimension: ");
    
    save(strcat(file_loc, file, file_type), 'cdims', '-append'); 
end

try
    load(strcat(file_loc, file, file_type), 'nn');
catch
    disp("Failed to load nn")
end

[m,n,p,~] = size(nn);

dx = m / cdims(1);
dy = n / cdims(2);
dz = p / cdims(3);

% Although the analytic computation of the uniform helical state yields 0,
% the derivatives are computed numerically so instead we shift energy with a zeropoint
% In this case, zeropoint_ee = zeropoint_eek22
if use_zeropoint
    zeropoint_ee = fFreeHelicalZeroPoint(m, n, p, k22, cdims(3), dz, chirality);
else
    zeropoint_ee = zeros(m-4,n-4,p-4);
end

if use_k22_only
    eesel = fFreek22(nn, k22, dx, dy, dz, chirality);
else
    eesel = fFree(nn, k11, k22, k33, k24, dx, dy, dz, chirality);
end

ee = eesel - zeropoint_ee;

sum(ee(:))

clear zeropoint_ee

nn = nn(3:end-2,3:end-2,3:end-2,:);

[m,n,p, ~] = size(nn);

p0 = 1;%floor(p/3);
p1 = p;%floor(p*2/3);

if xs_xy
    nn = nn(:,:,p0:p1,:);
    ee = ee(:,:,p0:p1);
    [m,n,p, ~] = size(nn);
end

disp('Computed ee') 

% Get total energy

total_energy = sum(sum(sum(ee))); % three times to get scalar

emin = min(min(min(ee)));
emax = max(max(max(ee)));

elevels = [];

if energy_iso
   elevels = [emin, emax]; 
end

cone = 0.015;

xy_midplane = ee(:,:,floor(p/2));
yz_midplane = ee(floor(m/2),:,:);
xz_midplane = ee(:,floor(n/2),:);

energy_per_voxel = total_energy / (m * n * p);
energy_per_pitch = total_energy / (cdims(1) * cdims(2) * cdims(3));
energy_per_voxel_per_pitch = energy_per_voxel / (cdims(1) * cdims(2) * cdims(3));

xx = linspace(1,m,m);
yy = linspace(1,n,n);
zz = linspace(1,p,p);
[xx,yy,zz] = ndgrid(xx,yy,zz);

for i=1:length(elevels)
    
    elevel = elevels(i);
    
    diffmag = abs(ee - elevel * ones(m,n,p));

    fv = isosurface(xx, yy, zz, diffmag, cone);
    disp('Created isosurface')

    thetad = 90;
    
    % Minimum blue, maximum red
    
    phid = 180 * abs(elevel-emax) / (emax-emin);

    theta = thetad * pi / 180;
    phi = phid * pi / 180;

    color = zeros(1,3);
    color(1) = phi/2/pi;
    color(2) = theta/(pi/2);
    color(3) = 2 - theta/(pi/2);
    
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

thetad = [0, 180];
phid = [0, 0];

epsilon = 0.22*ones(size(thetad));

theta = thetad/180*pi;
phi=phid/180*pi;

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

disp('interpolated')
% % normalize director field after interpolation
modnn = sqrt(nn_int(:,:,:,1).^2+nn_int(:,:,:,2).^2+nn_int(:,:,:,3).^2);
nn_int(:,:,:,1) = nn_int(:,:,:,1)./modnn;
nn_int(:,:,:,2) = nn_int(:,:,:,2)./modnn;
nn_int(:,:,:,3) = nn_int(:,:,:,3)./modnn;

disp('finished normalization')

nn2 = nn_int;
[ma,na,pa,~] = size(nn2);

for cone=1:length(theta)

    phixs = sin(theta(cone)).*cos(phi(cone));
    py = sin(theta(cone)).*sin(phi(cone));
    pz = cos(theta(cone));

    diffmag = sqrt((nn2(:,:,:,1)-phixs).^2+...
        (nn2(:,:,:,2)-py).^2+(nn2(:,:,:,3)-pz).^2);

    disp('Finished diffmag')

    fv = isosurface(xx,yy,zz,diffmag,epsilon(cone));

    disp('Created a isosurface')


    if numel(fv.vertices)>3

    % % old color scheme
    %     color = map(end-round(abs(theta(cone))/pi*(2^k-1)),:);

    % % HSV color scheme
    % %     linear
        color = zeros(1,3);
        color(1) = phi(cone)/2/pi;
        color(2) = theta(cone)/(pi/2);
        color(3) = 2 - theta(cone)/(pi/2);
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

    disp(['Made preimage (', cone, '/', length(theta), ')'])


end
    



if xs_xy
    
    x = linspace(1,m,m*ratexs);
    y = linspace(1,n,n*ratexs);
    [yy,xx] = meshgrid(y,x);
    ee_int2 = interp2(squeeze(xy_midplane),yy,xx,'cubic');
    
    cxy = colorbar;
    cxy.Label.String = "Free Energy";
    
    if use_zeropoint
        cxy.Label.String = strcat(cxy.Label.String, " (rel. bg)");
    end
    
    if use_k22_only
        cxy.Label.String = strcat("K22 ", cxy.Label.String);
    end
    S = surf(xx,yy,(p+1)/2*ones(size(xx)), ee_int2);
    %S = surf(xx, yy, ee_int2);
    S.EdgeColor = 'none';
    S.AmbientStrength = 0.3;
    S.FaceAlpha = xsection_alpha;
    
end

if xs_xz
   
    viewangle = [0,180];
    
    x = linspace(1,m,m*ratexs);
    z = linspace(1,p,p*ratexs);
    [zz,xx] = meshgrid(z,x);
    ee_int2 = interp2(squeeze(xz_midplane),zz,xx,'cubic');
    
    cxz = colorbar;
    cxz.Label.String = "Free Energy";
    
    if use_zeropoint
        cxz.Label.String = strcat(cxz.Label.String, " (rel. bg)");
    end
    
    if use_k22_only
        cxz.Label.String = strcat("K22 ", cxz.Label.String);
    end
    
    Sxz = surf(xx,(n+1)/2*ones(size(xx)),zz, ee_int2);
    Sxz.EdgeColor = 'none';
    %Sxz.AmbientStrength = 0.3;
    Sxz.FaceAlpha = xsection_alpha;
    
end

if xs_yz
    
    y = linspace(1,n,n*ratexs);
    z = linspace(1,p,p*ratexs);
    [zz,yy] = meshgrid(z,y);
    ee_int2 = interp2(squeeze(yz_midplane),zz,yy,'cubic');
    
    
    cyz = colorbar;
    cyz.Label.String = "Free Energy";
    
    if use_zeropoint
        cyz.Label.String = strcat(cyz.Label.String, " (rel. bg)");
    end
    
    if use_k22_only
        cyz.Label.String = strcat("K22 ", cyz.Label.String);
    end
    Syz = surf((m+1)/2*ones(size(yy)),yy,zz, ee_int2);
    Syz.EdgeColor = 'none';
    %Syz.AmbientStrength = 0.3;
    Syz.FaceAlpha = xsection_alpha;
    
end

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


function energy = fFree(nn, k11, k22, k33, k24, dx, dy, dz, chirality)

nx000 = nn(3:end-2,3:end-2,3:end-2,1);
ny000 = nn(3:end-2,3:end-2,3:end-2,2);
nz000 = nn(3:end-2,3:end-2,3:end-2,3);
nx100 = derivativeX(nn(:,:,:,1),dx);
nx010 = derivativeY(nn(:,:,:,1),dy);
nx001 = derivativeZ(nn(:,:,:,1),dz);
ny100 = derivativeX(nn(:,:,:,2),dx);
ny010 = derivativeY(nn(:,:,:,2),dy);
ny001 = derivativeZ(nn(:,:,:,2),dz);
nz100 = derivativeX(nn(:,:,:,3),dx);
nz010 = derivativeY(nn(:,:,:,3),dy);
nz001 = derivativeZ(nn(:,:,:,3),dz);


energy = (nx001.^2+nx010.^2+nx100.^2+ny001.^2+ny010.^2+ny100.^2+nz001.^2+(-1+...
k11).*(nx100+ny010+nz001).^2+nz010.^2+(-1+k33).*((nx000.*(-nx010+...
ny100)+nz000.*(ny001-nz010)).^2+(ny000.*(ny001-nz010)+nx000.*(nx001-...
nz100)).^2+(ny000.*(nx010-ny100)+nz000.*(nx001-nz100)).^2)+(-1+k22).*((-nx010+ny100).*nz000+...
nx000.*(-ny001+nz010)+ny000.*(nx001-nz100)).^2+...
nz100.^2-2.*(-1+k24).*(nx100.*ny010-nx010.*ny100+(nx100+ny010).*nz001-...
ny001.*nz010-nx001.*nz100)+4.*chirality.*k22.*pi.*((-nx010+ny100).*nz000+...
nx000.*(-ny001+nz010)+ny000.*(nx001-nz100)+pi))./2;

end

function energy = fFreek22(nn, k22, dx, dy, dz, chirality)

nx000 = nn(3:end-2,3:end-2,3:end-2,1);
ny000 = nn(3:end-2,3:end-2,3:end-2,2);
nz000 = nn(3:end-2,3:end-2,3:end-2,3);

nx010 = derivativeY(nn(:,:,:,1),dy);
nx001 = derivativeZ(nn(:,:,:,1),dz);
ny100 = derivativeX(nn(:,:,:,2),dx);

ny001 = derivativeZ(nn(:,:,:,2),dz);
nz100 = derivativeX(nn(:,:,:,3),dx);
nz010 = derivativeY(nn(:,:,:,3),dy);


energy = (k22.*(nx001.*ny000-nx010.*nz000+ny100.*nz000+nx000.*(-ny001+nz010)-...
ny000.*nz100+2.*chirality.*pi).^2)./2;

end

% n = (-sin(2pi z), cos(2pi z), 0), where z is dimensionless

function energy = fFreeHelicalZeroPoint(m, n, p, k22, cz, dz, chirality)

ztwist = linspace(0, 2 * pi * cz, p);
lateralx = linspace(1,1,m);
lateraly = linspace(1,1,n);
omega = ndgrid(lateralx, lateraly, ztwist);

nn = cat(4, -sin(omega), cos(omega), zeros(m,n,p));

nx000 = nn(3:end-2,3:end-2,3:end-2,1);
ny000 = nn(3:end-2,3:end-2,3:end-2,2);

nx001 = derivativeZ(nn(:,:,:,1),dz);
ny001 = derivativeZ(nn(:,:,:,2),dz);


energy = (k22.*(nx001.*ny000 - nx000.*ny001 + 2.*chirality.*pi).^2)./2;

end

function DxF = derivativeX(F, dx)
    % 4th order FD coefficients
    c_0 = 1/12/dx;
    c_1 = -2/3/dx;
    c_2 = 2/3/dx;
    c_3 = -1/12/dx;
    
    DxF = c_3 * F(5:end,3:end-2,3:end-2) + c_2 * F(4:end-1,3:end-2,3:end-2) +...
        c_1 * F(2:end-3,3:end-2,3:end-2) + c_0 * F(1:end-4,3:end-2,3:end-2);
end

function DyF = derivativeY(F, dy)
    % 4th order FD coefficients
    c_0 = 1/12/dy;
    c_1 = -2/3/dy;
    c_2 = 2/3/dy;
    c_3 = -1/12/dy;
    
    DyF = c_3 * F(3:end-2,5:end,3:end-2) + c_2 * F(3:end-2,4:end-1,3:end-2) +...
        c_1 * F(3:end-2,2:end-3,3:end-2) + c_0 * F(3:end-2,1:end-4,3:end-2);
end

function DzF = derivativeZ(F, dz)
    % 4th order FD coefficients
    c_0 = 1/12/dz;
    c_1 = -2/3/dz;
    c_2 = 2/3/dz;
    c_3 = -1/12/dz;
    
    DzF = c_3 * F(3:end-2,3:end-2,5:end) + c_2 * F(3:end-2,3:end-2,4:end-1) +...
        c_1 * F(3:end-2,3:end-2,2:end-3) + c_0 * F(3:end-2,3:end-2,1:end-4);
end
