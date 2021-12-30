% Compute the hopf index with default interpolation rate sample = 1 and interpolation
% method = 'cubic'
% -----------
% Modified by Darian to be more considerate to space and use 4th order FD
% -----------

function Q = hopf_index(nn, varargin)

Defaults = {1, 'cubic'};
idx = ~cellfun('isempty',varargin);
Defaults(idx) = varargin(idx);

sample = Defaults{1};
method = Defaults{2};

tic

size(nn);
[m,n,p,~] = size(nn);


% Original grid
xx = linspace(1,m,m);
yy = linspace(1,n,n);
zz = linspace(1,p,p);

% Refined grid to accomodate for central difference and
% interpolation
xxa = linspace(1,m,m*sample);
% xxa = [xxa(1)+xxa(1)-xxa(2) xxa xxa(end)+xxa(2)-xxa(1)];
yya = linspace(1,n,n*sample);
% yya = [yya(1)+yya(1)-yya(2) yya yya(end)+yya(2)-yya(1)];
zza = linspace(1,p,p*sample);
% zza = [zza(1)+zza(1)-zza(2) zza zza(end)+zza(2)-zza(1)];

[xx,yy,zz] = ndgrid(xx,yy,zz);
[xxa,yya,zza] = ndgrid(xxa,yya,zza);

nn_int(:,:,:,1) = interp3(yy,xx,zz,nn(:,:,:,1),yya,xxa,zza,method);
nn_int(:,:,:,2) = interp3(yy,xx,zz,nn(:,:,:,2),yya,xxa,zza,method);
nn_int(:,:,:,3) = interp3(yy,xx,zz,nn(:,:,:,3),yya,xxa,zza,method);
% clear xxa yya zza xx yy zz

% % normalize director field after interpolation
modnn = sqrt(nn_int(:,:,:,1).^2+nn_int(:,:,:,2).^2+nn_int(:,:,:,3).^2);
nn_int(:,:,:,1) = nn_int(:,:,:,1)./modnn;
nn_int(:,:,:,2) = nn_int(:,:,:,2)./modnn;
nn_int(:,:,:,3) = nn_int(:,:,:,3)./modnn;
clear modnn

% 4th order FD coefficients
c_0 = 1/12;
c_1 = -2/3;
c_2 = 2/3;
c_3 = -1/12;

% Creating B field Bi=E(ijk)nn.(D_j nn cross D_k nn)

partial_nn = zeros(m*sample,n*sample,p*sample,3,3); 
% coordinates, partial_i, n_j

partial_nn(3:end-2,:,:,1,:) = c_3 * nn_int(5:end,:,:,:) + c_2 * nn_int(4:end-1,:,:,:) + c_1 * nn_int(2:end-3,:,:,:) + c_0 * nn_int(1:end-4,:,:,:);
partial_nn(:,3:end-2,:,2,:) = c_3 * nn_int(:,5:end,:,:) + c_2 * nn_int(:,4:end-1,:,:) + c_1 * nn_int(:,2:end-3,:,:) + c_0 * nn_int(:,1:end-4,:,:);
partial_nn(:,:,3:end-2,3,:) = c_3 * nn_int(:,:,5:end,:) + c_2 * nn_int(:,:,4:end-1,:) + c_1 * nn_int(:,:,2:end-3,:) + c_0 * nn_int(:,:,1:end-4,:);


B = zeros(m*sample,n*sample,p*sample,3);

B(:,:,:,1) = 2*dot(nn_int,cross(squeeze(partial_nn(:,:,:,2,:)),squeeze(partial_nn(:,:,:,3,:)),4),4);
B(:,:,:,2) = 2*dot(nn_int,cross(squeeze(partial_nn(:,:,:,3,:)),squeeze(partial_nn(:,:,:,1,:)),4),4);
B(:,:,:,3) = 2*dot(nn_int,cross(squeeze(partial_nn(:,:,:,1,:)),squeeze(partial_nn(:,:,:,2,:)),4),4);
clear partial_nn nn_int
% Construct A. Choose a gauge such that Az = 0;
A = zeros(m*sample,n*sample,p*sample,3);

A(:,:,:,2) = -cumtrapz(B(:,:,:,1),3);
A(:,:,:,1) = cumtrapz(B(:,:,:,2),3);

% partialx Ay - partialy Ax = Bz
Bz = B(:,:,:,3);

% Calculate c1
pxAy = c_3*A(5:end,:,:,2)+c_2*A(4:end-1,:,:,2)+c_1*A(2:end-3,:,:,2)+ c_0*A(1:end-4,:,:,2);
pyAx = c_3*A(:,5:end,:,1)+c_2*A(:,4:end-1,:,1)+c_1*A(:,2:end-3,:,1)+ c_0*A(:,1:end-4,:,1);
pxAy = cat(1,zeros(2,n*sample,p*sample),pxAy,zeros(2,n*sample,p*sample));
pyAx = cat(2,zeros(m*sample,2,p*sample),pyAx,zeros(m*sample,2,p*sample));

c1 = cumtrapz(-Bz + pxAy - pyAx,2);
clear pxAy pyAx

% % Check B and curl_A difference and Calculate Hopf index

% % Without C1
% Hopf_index1 = trapz(trapz(trapz(dot(B,A,4))))/(8*pi)^2

% % With C1
A(:,:,:,1) = A(:,:,:,1) + c1;

Hopf_index2 = trapz(trapz(trapz(dot(B,A,4))))/(8*pi)^2;

toc

Q = Hopf_index2;

fprintf('Hopf index = %0.20e.',Q);

%hold all
%Qdensity = dot(B,A,4)/(8*pi)^2;
%Qdensity(xxa > (m+1)/2) = nan;
%iso = [0.1 0.5 1 3 6]*1E-5;
%for i = 1:length(iso)
%fv = isosurface(xxa,yya,zza,Qdensity,iso(i));  
%patch(fv,...
%        'FaceColor',[iso(i)/max(iso) 0 0],...
%        'FaceLighting','gouraud',...
%        'FaceAlpha',0.5,...
%        'EdgeColor','none',...
%        'SpecularStrength',0.7,...
%        'AmbientStrength',0.5);
%end
%camlight
%axis equal
%xlim([1 m])
%ylim([1 n])
%zlim([1 p])

end