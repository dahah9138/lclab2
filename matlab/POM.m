% Code originally written by Benny
% Modified by Darian to be more efficient with memory Mar. 16, 2021

clear
clc
% clf
figure

file = 'voltage_dimer_physical.mat';
%file = 'Output/toron.mat';
crossed = 1;
blur = 1;
scan_ratio = 1; % Ratio of the cell to include in POM

pitch = 4; % pitch in microns
cell_dop = 5;
dop = 5; % thickness over pitch of simulated heliknoton
layers_to_add_below = cell_dop - dop;
layers_to_add_above =0;

try
    load(file,'nn');
catch
   
    % See what a uniform radial structure looks like
    sample = 2;
    
    m = sample * 40;
    n = sample * 40;
    p = sample * 10;
    
    x = linspace(-1,1, m);
    y = linspace(-1,1, n);
    z = linspace(-0.1, 0.1, p);
    
    [x,y,z] = ndgrid(x,y,z);
    
    [~,theta,~] = cart2pol(x,y,z);
    
    nn = zeros([size(theta), 3]);
    
    nn(:,:,:,1) = cos(theta);
    nn(:,:,:,2) = sin(theta);
    nn(:,:,:,3) = 0;
    
    clear theta x y z
    
end

% % Upsample and interpolate
sampleXY = 4;
sampleZ = 1;

% % POM setup
waveplate = 0; % w/ or w/out a 530 nm wave retardation plate; https://www.olympus-lifescience.com/en/microscope-resource/primer/techniques/polarized/firstorderplate/
th = crossed * 90*pi/180; % angle between polarizer and analyzer

% % cell thickness 
thickness = pitch*(dop+layers_to_add_above+layers_to_add_below); % cell thickness in microns
lambda=[650,550,450]*10^-9; % RGB wavelength
intensity=[1 0.6 0.2]*1;  % intensity factor accounting for lamp spectrum and sensitivity anisotropy at different wavelengths
gamma = 1;
%lambda=[650,510,475]*10^-9;

% % Material birefringence
% % 5CB
n0 = 1.58;
ne = 1.77;
% E7
% n0 = 1.52;
% ne = 1.73;
% AMLC
% n0=1.47;
% ne=1.55;


% % Interpolation
[m,n,p,d] = size(nn);
mp = m*sampleXY;
np = n*sampleXY;
pp = p*sampleZ;

nodes_in_layer = cast(pp / dop, 'uint16');

xx = linspace(1,m,mp);
yy = linspace(1,n,np);
zz = linspace(1,p,pp);
[yy,xx,zz] = meshgrid(yy,xx,zz);
 
nn_int(:,:,:,1) = interp3(nn(:,:,:,1),yy,xx,zz,'cubic');
nn_int(:,:,:,2) = interp3(nn(:,:,:,2),yy,xx,zz,'cubic');
nn_int(:,:,:,3) = interp3(nn(:,:,:,3),yy,xx,zz,'cubic');

% % normalize director field after interpolation
modnn = sqrt(nn_int(:,:,:,1).^2+nn_int(:,:,:,2).^2+nn_int(:,:,:,3).^2);
nn_int(:,:,:,1) = nn_int(:,:,:,1)./modnn;
nn_int(:,:,:,2) = nn_int(:,:,:,2)./modnn;
nn_int(:,:,:,3) = nn_int(:,:,:,3)./modnn;

clear nn
nn = nn_int;
dz=thickness*10^(-6)/pp;
phi = atan2(nn(:,:,:,2),nn(:,:,:,1));
theta = pi/2 - atan2(nn(:,:,:,3),sqrt(nn(:,:,:,1).^2 + nn(:,:,:,2).^2));

C = zeros(mp,np,3);

tic

% Eo_0 = [1,0]^t for each node in xy plane for each color
Eo=repmat([1;0],1,mp,np,3);

% Jones matrix for each (x,y,~) stack
M = zeros(2,2,mp,np);


nodes_in_layers = nodes_in_layer*layers_to_add_below;
dPhi = 1 / cast(nodes_in_layer - 1, 'double');

% Iterate through rgb colors
for l = 1:3
    
    % Iterate through a single planar twisting node in xy plane
    % result is [E0(1,ell); E0(2,ell)] = Mell(p)Mell(p-1)...Mell(1)[1;0]
    for t = 1:(nodes_in_layers-1)
        
        delta0=2*pi/lambda(l)*dz*n0;
        deltaE=2*pi./lambda(l)*dz*ne;
        ang = cast(2 * pi * dPhi * cast(t - 1, 'double') + pi/2, 'double');
        
        M(1,1,:,:)=cos(ang)^2*exp(1i*deltaE)+sin(ang)^2*exp(1i*delta0); 
        M(1,2,:,:)=sin(ang)*cos(ang)*(exp(1i*deltaE)-exp(1i*delta0));
        M(2,1,:,:)=M(1,2,:,:);
        M(2,2,:,:)=sin(ang)^2*exp(1i*deltaE)+cos(ang)^2*exp(1i*delta0);
        
        E1=squeeze(Eo(1,:,:,l));
        E2=squeeze(Eo(2,:,:,l));
        Eo(1,:,:,l)=squeeze(M(1,1,:,:)).*E1+squeeze(M(1,2,:,:)).*E2;
        Eo(2,:,:,l)=squeeze(M(2,1,:,:)).*E1+squeeze(M(2,2,:,:)).*E2;
    end

end


% Iterate through rgb colors
for l = 1:3
    
    % Iterate through xy cross sections multiplying them together
    % result is [E0(1,ell); E0(2,ell)] = Mell(p)Mell(p-1)...Mell(1)[1;0]
    for t = 1:int32(pp*scan_ratio)


        delta0=2*pi/lambda(l)*dz*n0;
        ne_th=ne*n0./sqrt(n0^2*(sin(theta(:,:,t))).^2+ne^2*(cos(theta(:,:,t))).^2);
        deltaE=2*pi./lambda(l)*dz*ne_th;

        M(1,1,:,:)=cos(phi(:,:,t)).^2.*exp(1i*deltaE)+sin(phi(:,:,t)).^2*exp(1i*delta0); 
        M(1,2,:,:)=sin(phi(:,:,t)).*cos(phi(:,:,t)).*(exp(1i*deltaE)-exp(1i*delta0));
        M(2,1,:,:)=sin(phi(:,:,t)).*cos(phi(:,:,t)).*(exp(1i*deltaE)-exp(1i*delta0)); 
        M(2,2,:,:)=sin(phi(:,:,t)).^2.*exp(1i*deltaE)+cos(phi(:,:,t)).^2*exp(1i*delta0);

        E1=squeeze(Eo(1,:,:,l));
        E2=squeeze(Eo(2,:,:,l));
        Eo(1,:,:,l)=squeeze(M(1,1,:,:)).*E1+squeeze(M(1,2,:,:)).*E2;
        Eo(2,:,:,l)=squeeze(M(2,1,:,:)).*E1+squeeze(M(2,2,:,:)).*E2;

    end

end


nodes_in_layers = nodes_in_layer*layers_to_add_above;
dPhi = 1 / cast(nodes_in_layer - 1, 'double');

% Iterate through rgb colors
for l = 1:3
    
    % Iterate through a single planar twisting node in xy plane
    % result is [E0(1,ell); E0(2,ell)] = Mell(p)Mell(p-1)...Mell(1)[1;0]
    for t = 1:(nodes_in_layers-1)
        
        delta0=2*pi/lambda(l)*dz*n0;
        deltaE=2*pi./lambda(l)*dz*ne;
        ang = cast(2 * pi * dPhi * cast(t, 'double') + pi/2, 'double');
        
        M(1,1,:,:)=cos(ang)^2*exp(1i*deltaE)+sin(ang)^2*exp(1i*delta0); 
        M(1,2,:,:)=sin(ang)*cos(ang)*(exp(1i*deltaE)-exp(1i*delta0));
        M(2,1,:,:)=M(1,2,:,:);
        M(2,2,:,:)=sin(ang)^2*exp(1i*deltaE)+cos(ang)^2*exp(1i*delta0);

        E1=squeeze(Eo(1,:,:,l));
        E2=squeeze(Eo(2,:,:,l));
        Eo(1,:,:,l)=squeeze(M(1,1,:,:)).*E1+squeeze(M(1,2,:,:)).*E2;
        Eo(2,:,:,l)=squeeze(M(2,1,:,:)).*E1+squeeze(M(2,2,:,:)).*E2;

    end

end

clear theta phi

switch waveplate
    case 1
     % 530nm full waveplate
    for l = [1 3]

    thick=5.8889e-05*7;
    de_wp=2*pi/lambda(l)*thick*1.55338;
    do_ep=2*pi/lambda(l)*thick*1.54425;

    m(1,1)=0.5*(exp(1i*de_wp)+exp(1i*do_ep)); 
    m(1,2)=0.5*(exp(1i*de_wp)-exp(1i*do_ep)); 
    m(2,1)=0.5*(exp(1i*de_wp)-exp(1i*do_ep));  
    m(2,2)=0.5*(exp(1i*de_wp)+exp(1i*do_ep)); 

    E1=squeeze(Eo(1,:,:,l));
    E2=squeeze(Eo(2,:,:,l));
    Eo(1,:,:,l)=m(1,1).*E1+m(1,2).*E2;
    Eo(2,:,:,l)=m(2,1).*E1+m(2,2).*E2;
    end
    case 0
end

for l=1:3 
    
    I(:,:)=abs(Eo(1,:,:,l)*cos(th)+Eo(2,:,:,l)*sin(th)).^2;
    I = imgaussfilt(I,blur*sampleXY);
    C(:,:,l)=intensity(l)*I;
    C(:,:,l)=C(:,:,l).^gamma;
end

toc

% f=figure;
clf
imagesc(permute(C,[2 1 3]))
axis equal off
xlabel('x')
ylabel('y')
set(gca,'YDir','normal')

% switch waveplate
%     case 1
%         title(['POM 530wp, sampleXY=' num2str(sampleXY) ', sampleZ=' num2str(sampleZ)])
%     case 0
%         title(['POM, sampleXY=' num2str(sampleXY) ', sampleZ=' num2str(sampleZ)])
% end
% 
% switch waveplate
%     case 1
%         saveas(f,['POM 530wp, sampleXY' num2str(sampleXY) ', sampleZ' num2str(sampleZ)],'bmp')
%     case 0
%         saveas(f,['POM, sampleXY' num2str(sampleXY) ', sampleZ' num2str(sampleZ)],'bmp')
% end


