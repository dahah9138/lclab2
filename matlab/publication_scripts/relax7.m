%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RELAX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical User Interface to relax Liquid Crystal (LC) director

function relax7
% [relax] Brief description:
%   Run relax to generate the Graphical User Interface (GUI).
%   
% setappdata(hgui,'fieldname',value) creates or replaces a structure field
% data = getappdata(hgui) retrieves the data structure
% fieldname = getappdata(hgui,'fieldname') retrieves a specific field
%   

% enforces singleton (only one GUI at a time)
hgui = findall(0, 'Name', mfilename);
if isempty(hgui) % if no GUI existing, run make_fig_fcn
    % Create and hide the GUI figure as it is being constructed.
    hgui = make_fig_fcn;
end

% raise GUI and make visible to user
figure(hgui)
hgui.Visible = 'on';

%% %%%%%%%%%%%%%%%%%%%% Utility functions for RELAX %%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The physics!

function frelax(hgui,report,num_it,num_rep)
% [frelax] Brief description:
%   Shell for relaxation routien. 

% Get the data structure from the appdata container
data = getappdata(hgui);

% populate variable to be used with information from the data structure
nn = data.nn;                       % director data (4D array)
vv = data.vv;                       % voltage profile data (3D array)
m = str2double(data.edit_nx);       % grid points
n = str2double(data.edit_ny);
p = str2double(data.edit_nz);
mnp = m*n*p;
dimx = str2double(data.edit_dimx);  % size of cell in um
dimy = str2double(data.edit_dimy);
dimz = str2double(data.edit_dimz);
dx = dimx/m;                        % step size in um
dy = dimy/n;
dz = dimz/p;
k11 = str2double(data.edit_k11);    % elastic constants
k22 = str2double(data.edit_k22);
k33 = str2double(data.edit_k33);
k24 = str2double(data.edit_k24);
pitch = str2double(data.edit_pitch); % pitch value
q0 = 2*pi/pitch;                     % chiral wave number
eper = str2double(data.edit_epsper); % dielectric constants
ea = str2double(data.edit_epspar)-eper;
e0 = 8.854e-12;                      % permittivity of free space
hx = str2double(data.edit_mx);       % applied magnetic fields
hy = str2double(data.edit_my);
hz = str2double(data.edit_mz);
u0 = 1.256e-6;                       % permeability of free space

it = str2double(data.edit_iterations); % number of iterations to relax
avefg = zeros(it,1); % initialize the functional derivative array

% gain parameter for relaxation time step
gain = str2double(data.edit_gain);
% MaxStableTimeStep derived from heat equation
MSTS = min([dx dy dz])^2/2/max([k11 k22 k33 k24])*gain;

% Boundry condition vector: ex. [0 0 0] -> no periodic boundarys
BCs = [data.rb_periodicBCsTB data.rb_periodicBCsLR data.rb_periodicBCsFB];

% Get start time
t1 = tic;
for tt = 1:it
    % apply update formula to all lattice points
    % calcualte the functional derivatives (elastic[fsi] and electric[fei])
    [fsx, fsy, fsz] = fFsnn3D(nn,dx,dy,dz,k11,k22,k33,k24,q0);
    [fex, fey, fez] = fFenn3D(nn,vv,dx,dy,dz,ea,e0);
    [fmx, fmy, fmz] = fFmnn3D(m,n,p,hx,hy,hz,u0);
    fg = (fsx+fsy+fsz)-(fex+fey+fez)+(fmx+fmy+fmz);
    avefg(tt,1) = sum(fg(:))/mnp;
    % calculate new director for all lattice points
    nn(2:end-1,2:end-1,2:end-1,1) = nn(2:end-1,2:end-1,2:end-1,1)-...
        MSTS/2*(fsx-fex+fmx);
    nn(2:end-1,2:end-1,2:end-1,2) = nn(2:end-1,2:end-1,2:end-1,2)-...
        MSTS/2*(fsy-fey+fmy);
    nn(2:end-1,2:end-1,2:end-1,3) = nn(2:end-1,2:end-1,2:end-1,3)-...
        MSTS/2*(fsz-fez+fmz);
    % fix director here
    
% % % for relaxing only the interior cylindrical region
% % [xx,yy,~] = meshgrid(1:n,1:m,1:p);
% % xx = xx-(m+1)/2;
% % yy = yy-(n+1)/2;
% % nx = nn(:,:,:,1);
% % ny = nn(:,:,:,2);
% % nz = nn(:,:,:,3);
% % nx(sqrt(xx.^2+yy.^2)>(m+1)/2-2) = 0;
% % ny(sqrt(xx.^2+yy.^2)>(m+1)/2-2) = 0;
% % nz(sqrt(xx.^2+yy.^2)>(m+1)/2-2) = -1; % this can be 1 or -1
% % nn = cat(4,nx,ny,nz);
    
    % Normalize director
    mag = repmat(sqrt(sum((nn).^2,4)),1,1,1,3);
    nn = nn./mag;
    
    % BCs
    switch num2str(BCs)
        case '0  0  0'
        case '0  0  1'% f/b walls Periodic
            nn(:,1,:,:) = nn(:,n-1,:,:);nn(:,n,:,:) = nn(:,2,:,:); 
            vv(:,2,:,:) = vv(:,end-2,:,:);vv(:,end-1,:,:) = vv(:,3,:,:); 
            vv(:,1,:,:) = vv(:,end-1,:,:);vv(:,end,:,:) = vv(:,2,:,:); 
        case '0  1  0'%l/r walls Periodic
            nn(1,:,:,:) = nn(m-1,:,:,:);nn(m,:,:,:) = nn(2,:,:,:); 
            vv(2,:,:,:) = vv(end-2,:,:,:);vv(end-1,:,:,:) = vv(3,:,:,:); 
            vv(1,:,:,:) = vv(end-1,:,:,:);vv(end,:,:,:) = vv(2,:,:,:); 
        case '0  1  1'%f/b walls Periodic %l/r walls Periodic
            nn(:,1,:,:) = nn(:,n-1,:,:);nn(:,n,:,:) = nn(:,2,:,:); 
            nn(1,:,:,:) = nn(m-1,:,:,:);nn(m,:,:,:) = nn(2,:,:,:); 
            vv(:,2,:,:) = vv(:,end-2,:,:);vv(:,end-1,:,:) = vv(:,3,:,:); 
            vv(2,:,:,:) = vv(end-2,:,:,:);vv(end-1,:,:,:) = vv(3,:,:,:); 
            vv(:,1,:,:) = vv(:,end-1,:,:);vv(:,end,:,:) = vv(:,2,:,:); 
            vv(1,:,:,:) = vv(end-1,:,:,:);vv(end,:,:,:) = vv(2,:,:,:); 
        case '1  0  0'%t/b Periodic
            nn(:,:,1,:) = nn(:,:,p-1,:);nn(:,:,p,:) = nn(:,:,2,:); 
            vv(:,:,2,:) = vv(:,:,end-2,:);vv(:,:,end-1,:) = vv(:,:,3,:); 
            vv(:,:,1,:) = vv(:,:,end-1,:);vv(:,:,end,:) = vv(:,:,2,:); 
        case '1  0  1'%f/b walls Periodic %t/b Periodic
            nn(:,1,:,:) = nn(:,n-1,:,:);nn(:,n,:,:) = nn(:,2,:,:); 
            nn(:,:,1,:) = nn(:,:,p-1,:);nn(:,:,p,:) = nn(:,:,2,:); 
            vv(:,2,:,:) = vv(:,end-2,:,:);vv(:,end-1,:,:) = vv(:,3,:,:); 
            vv(:,:,2,:) = vv(:,:,end-2,:);vv(:,:,end-1,:) = vv(:,:,3,:); 
            vv(:,1,:,:) = vv(:,end-1,:,:);vv(:,end,:,:) = vv(:,2,:,:); 
            vv(:,:,1,:) = vv(:,:,end-1,:);vv(:,:,end,:) = vv(:,:,2,:); 
        case '1  1  0'%l/r walls Periodic %t/b Periodic
            nn(1,:,:,:) = nn(m-1,:,:,:);nn(m,:,:,:) = nn(2,:,:,:); 
            nn(:,:,1,:) = nn(:,:,p-1,:);nn(:,:,p,:) = nn(:,:,2,:); 
            vv(2,:,:,:) = vv(end-2,:,:,:);vv(end-1,:,:,:) = vv(3,:,:,:); 
            vv(:,:,2,:) = vv(:,:,end-2,:);vv(:,:,end-1,:) = vv(:,:,3,:); 
            vv(1,:,:,:) = vv(end-1,:,:,:);vv(end,:,:,:) = vv(2,:,:,:); 
            vv(:,:,1,:) = vv(:,:,end-1,:);vv(:,:,end,:) = vv(:,:,2,:); 
        case '1  1  1'%f/b walls Periodic %l/r walls Periodic %t/b Periodic
            nn(:,1,:,:) = nn(:,n-1,:,:);nn(:,n,:,:) = nn(:,2,:,:); 
            nn(1,:,:,:) = nn(m-1,:,:,:);nn(m,:,:,:) = nn(2,:,:,:); 
            nn(:,:,1,:) = nn(:,:,p-1,:);nn(:,:,p,:) = nn(:,:,2,:); 
            vv(:,2,:,:) = vv(:,end-2,:,:);vv(:,end-1,:,:) = vv(:,3,:,:); 
            vv(2,:,:,:) = vv(end-2,:,:,:);vv(end-1,:,:,:) = vv(3,:,:,:); 
            vv(:,:,2,:) = vv(:,:,end-2,:);vv(:,:,end-1,:) = vv(:,:,3,:); 
            vv(:,1,:,:) = vv(:,end-1,:,:);vv(:,end,:,:) = vv(:,2,:,:); 
            vv(1,:,:,:) = vv(end-1,:,:,:);vv(end,:,:,:) = vv(2,:,:,:); 
            vv(:,:,1,:) = vv(:,:,end-1,:);vv(:,:,end,:) = vv(:,:,2,:); 
    end % switch num2str(BCs)
    
    % update v for all lattice points
    vv = fvoltageUpdate(nn,vv,dx,dy,dz,ea,eper);
    setappdata(hgui,'nn',nn)
    setappdata(hgui,'vv',vv)
    
    % update iterations
    t2 = toc(t1);
    if(t2 > 0.05)
      htext_it = findobj('tag','text_it');
      set(htext_it,'String',tt + (report-1)*num_it/num_rep)
      drawnow
      t1 = tic;
    end % if(t2 > 0.05)
end % for tt = 1:it

% Calculate Ms (change in average fg over timestep and total time)
if data.rb_ms
    delta = [nan; avefg(1:end-1)-avefg(2:end)];
    Ms = [delta/MSTS linspace(0,MSTS*it,it)'];
    haxis_ms = findobj('tag','axis_ms');
    haxis = findobj('tag','axis');
    set(hgui, 'currentaxes', haxis_ms)
    plot(Ms(:,2)/MSTS,Ms(:,1),'k.-');axis tight;grid on
    set(haxis_ms,'tag','axis_ms')
    set(hgui, 'currentaxes', haxis)
end % if data.rb_ms

function [fsx, fsy, fsz] = fFsnn3D(nn,dx,dy,dz,k11,k22,k33,k24,q0)
%fFsnn3D calculates the functional derivatives (used in frelax3D)
nx000 = nn(2:end-1,2:end-1,2:end-1,1);
ny000 = nn(2:end-1,2:end-1,2:end-1,2);
nz000 = nn(2:end-1,2:end-1,2:end-1,3);
nx100 = (nn(3:end,2:end-1,2:end-1,1)-nn(1:end-2,2:end-1,2:end-1,1))/(2*dx);
nx010 = (nn(2:end-1,3:end,2:end-1,1)-nn(2:end-1,1:end-2,2:end-1,1))/(2*dy);
nx001 = (nn(2:end-1,2:end-1,3:end,1)-nn(2:end-1,2:end-1,1:end-2,1))/(2*dz);
ny100 = (nn(3:end,2:end-1,2:end-1,2)-nn(1:end-2,2:end-1,2:end-1,2))/(2*dx);
ny010 = (nn(2:end-1,3:end,2:end-1,2)-nn(2:end-1,1:end-2,2:end-1,2))/(2*dy);
ny001 = (nn(2:end-1,2:end-1,3:end,2)-nn(2:end-1,2:end-1,1:end-2,2))/(2*dz);
nz100 = (nn(3:end,2:end-1,2:end-1,3)-nn(1:end-2,2:end-1,2:end-1,3))/(2*dx);
nz010 = (nn(2:end-1,3:end,2:end-1,3)-nn(2:end-1,1:end-2,2:end-1,3))/(2*dy);
nz001 = (nn(2:end-1,2:end-1,3:end,3)-nn(2:end-1,2:end-1,1:end-2,3))/(2*dz);
nx200 = (nn(3:end,2:end-1,2:end-1,1)+nn(1:end-2,2:end-1,2:end-1,1)-2*nx000)/dx^2;
nx020 = (nn(2:end-1,3:end,2:end-1,1)+nn(2:end-1,1:end-2,2:end-1,1)-2*nx000)/dy^2;
nx002 = (nn(2:end-1,2:end-1,3:end,1)+nn(2:end-1,2:end-1,1:end-2,1)-2*nx000)/dz^2;
ny200 = (nn(3:end,2:end-1,2:end-1,2)+nn(1:end-2,2:end-1,2:end-1,2)-2*ny000)/dx^2;
ny020 = (nn(2:end-1,3:end,2:end-1,2)+nn(2:end-1,1:end-2,2:end-1,2)-2*ny000)/dy^2;
ny002 = (nn(2:end-1,2:end-1,3:end,2)+nn(2:end-1,2:end-1,1:end-2,2)-2*ny000)/dz^2;
nz200 = (nn(3:end,2:end-1,2:end-1,3)+nn(1:end-2,2:end-1,2:end-1,3)-2*nz000)/dx^2;
nz020 = (nn(2:end-1,3:end,2:end-1,3)+nn(2:end-1,1:end-2,2:end-1,3)-2*nz000)/dy^2;
nz002 = (nn(2:end-1,2:end-1,3:end,3)+nn(2:end-1,2:end-1,1:end-2,3)-2*nz000)/dz^2;
nx110 = (nn(3:end,3:end,2:end-1,1)-nn(3:end,1:end-2,2:end-1,1)-nn(1:end-2,3:end,2:end-1,1)+nn(1:end-2,1:end-2,2:end-1,1))/(4*dx*dy);
ny110 = (nn(3:end,3:end,2:end-1,2)-nn(3:end,1:end-2,2:end-1,2)-nn(1:end-2,3:end,2:end-1,2)+nn(1:end-2,1:end-2,2:end-1,2))/(4*dx*dy);
nz110 = (nn(3:end,3:end,2:end-1,3)-nn(3:end,1:end-2,2:end-1,3)-nn(1:end-2,3:end,2:end-1,3)+nn(1:end-2,1:end-2,2:end-1,3))/(4*dx*dy);
nx101 = (nn(3:end,2:end-1,3:end,1)-nn(3:end,2:end-1,1:end-2,1)-nn(1:end-2,2:end-1,3:end,1)+nn(1:end-2,2:end-1,1:end-2,1))/(4*dx*dz);
ny101 = (nn(3:end,2:end-1,3:end,2)-nn(3:end,2:end-1,1:end-2,2)-nn(1:end-2,2:end-1,3:end,2)+nn(1:end-2,2:end-1,1:end-2,2))/(4*dx*dz);
nz101 = (nn(3:end,2:end-1,3:end,3)-nn(3:end,2:end-1,1:end-2,3)-nn(1:end-2,2:end-1,3:end,3)+nn(1:end-2,2:end-1,1:end-2,3))/(4*dx*dz);
nx011 = (nn(2:end-1,3:end,3:end,1)-nn(2:end-1,3:end,1:end-2,1)-nn(2:end-1,1:end-2,3:end,1)+nn(2:end-1,1:end-2,1:end-2,1))/(4*dy*dz);
ny011 = (nn(2:end-1,3:end,3:end,2)-nn(2:end-1,3:end,1:end-2,2)-nn(2:end-1,1:end-2,3:end,2)+nn(2:end-1,1:end-2,1:end-2,2))/(4*dy*dz);
nz011 = (nn(2:end-1,3:end,3:end,3)-nn(2:end-1,3:end,1:end-2,3)-nn(2:end-1,1:end-2,3:end,3)+nn(2:end-1,1:end-2,1:end-2,3))/(4*dy*dz);

fsx = -(k33.*nx000.*nx001.^2)+(k24.*nx002)/2-k33.*nx000.^2.*nx002-k33.*nx000.*nx010.^2+(k24.*nx020)/2-k33.*nx000.^2.*nx020-k11.*nx200+(k24.*nx200)/2-k22.*nx002.*ny000.^2-k33.*nx020.*ny000.^2-2.*k22.*nx001.*ny000.*ny001+2.*k22.*nx000.*ny001.^2-k33.*nx000.*ny001.^2+k22.*nx000.*ny000.*ny002-k33.*nx000.*ny000.*ny002-2.*k33.*nx010.*ny000.*ny010+2.*k33.*ny000.*ny010.*ny100+k33.*nx000.*ny100.^2-k11.*ny110+k33.*nx000.^2.*ny110+k33.*ny000.^2.*ny110+2.*k22.*nx011.*ny000.*nz000-2.*k33.*nx011.*ny000.*nz000+k22.*nx010.*ny001.*nz000-k33.*nx010.*ny001.*nz000+k22.*nx001.*ny010.*nz000-k33.*nx001.*ny010.*nz000-k22.*nx000.*ny011.*nz000+k33.*nx000.*ny011.*nz000-2.*k22.*ny001.*ny100.*nz000+2.*k33.*ny001.*ny100.*nz000-k22.*ny000.*ny101.*nz000+k33.*ny000.*ny101.*nz000-k33.*nx002.*nz000.^2-k22.*nx020.*nz000.^2+k22.*ny110.*nz000.^2+k22.*nx010.*ny000.*nz001-k33.*nx010.*ny000.*nz001-k22.*ny000.*ny100.*nz001+k33.*ny000.*ny100.*nz001-2.*k33.*nx001.*nz000.*nz001+k22.*nx001.*ny000.*nz010-k33.*nx001.*ny000.*nz010-4.*k22.*nx000.*ny001.*nz010+2.*k33.*nx000.*ny001.*nz010-2.*k22.*nx010.*nz000.*nz010+3.*k22.*ny100.*nz000.*nz010-k33.*ny100.*nz000.*nz010+2.*k22.*nx000.*nz010.^2-k33.*nx000.*nz010.^2-k22.*nx000.*ny000.*nz011+k33.*nx000.*ny000.*nz011+k22.*nx000.*nz000.*nz020-k33.*nx000.*nz000.*nz020+3.*k22.*ny000.*ny001.*nz100-k33.*ny000.*ny001.*nz100-k22.*ny010.*nz000.*nz100+k33.*ny010.*nz000.*nz100+2.*k33.*nz000.*nz001.*nz100-2.*k22.*ny000.*nz010.*nz100+2.*k33.*ny000.*nz010.*nz100+k33.*nx000.*nz100.^2-k11.*nz101+k33.*nx000.^2.*nz101+k22.*ny000.^2.*nz101+k33.*nz000.^2.*nz101-k22.*ny000.*nz000.*nz110+k33.*ny000.*nz000.*nz110-2.*k22.*ny001.*q0+2.*k22.*nz010.*q0;
fsy = 2.*k33.*nx000.*nx010.*nx100-k11.*nx110+k33.*nx000.^2.*nx110+2.*k22.*nx001.^2.*ny000-k33.*nx001.^2.*ny000+k22.*nx000.*nx002.*ny000-k33.*nx000.*nx002.*ny000+k33.*nx010.^2.*ny000+k33.*nx110.*ny000.^2-2.*k22.*nx000.*nx001.*ny001-k33.*ny000.*ny001.^2+(k24.*ny002)/2-k22.*nx000.^2.*ny002-k33.*ny000.^2.*ny002-k11.*ny020+(k24.*ny020)/2-2.*k33.*nx000.*nx100.*ny100-k33.*ny000.*ny100.^2+(k24.*ny200)/2-k33.*nx000.^2.*ny200-k33.*ny000.^2.*ny200-2.*k22.*nx001.*nx010.*nz000+2.*k33.*nx001.*nx010.*nz000-k22.*nx000.*nx011.*nz000+k33.*nx000.*nx011.*nz000-k22.*nx101.*ny000.*nz000+k33.*nx101.*ny000.*nz000+k22.*nx100.*ny001.*nz000-k33.*nx100.*ny001.*nz000+k22.*nx001.*ny100.*nz000-k33.*nx001.*ny100.*nz000+2.*k22.*nx000.*ny101.*nz000-2.*k33.*nx000.*ny101.*nz000+k22.*nx110.*nz000.^2-k33.*ny002.*nz000.^2-k22.*ny200.*nz000.^2-k22.*nx000.*nx010.*nz001+k33.*nx000.*nx010.*nz001+k22.*nx000.*ny100.*nz001-k33.*nx000.*ny100.*nz001-2.*k33.*ny001.*nz000.*nz001+3.*k22.*nx000.*nx001.*nz010-k33.*nx000.*nx001.*nz010-k22.*nx100.*nz000.*nz010+k33.*nx100.*nz000.*nz010+2.*k33.*nz000.*nz001.*nz010+k33.*ny000.*nz010.^2-k11.*nz011+k22.*nx000.^2.*nz011+k33.*ny000.^2.*nz011+k33.*nz000.^2.*nz011-4.*k22.*nx001.*ny000.*nz100+2.*k33.*nx001.*ny000.*nz100+k22.*nx000.*ny001.*nz100-k33.*nx000.*ny001.*nz100+3.*k22.*nx010.*nz000.*nz100-k33.*nx010.*nz000.*nz100-2.*k22.*ny100.*nz000.*nz100-2.*k22.*nx000.*nz010.*nz100+2.*k33.*nx000.*nz010.*nz100+2.*k22.*ny000.*nz100.^2-k33.*ny000.*nz100.^2-k22.*nx000.*ny000.*nz101+k33.*nx000.*ny000.*nz101-k22.*nx000.*nz000.*nz110+k33.*nx000.*nz000.*nz110+k22.*ny000.*nz000.*nz200-k33.*ny000.*nz000.*nz200+2.*k22.*nx001.*q0-2.*k22.*nz100.*q0;
fsz = 2.*k33.*nx000.*nx001.*nx100-k11.*nx101+k33.*nx000.^2.*nx101-2.*k22.*nx001.*nx010.*ny000+2.*k33.*nx001.*nx010.*ny000-k22.*nx000.*nx011.*ny000+k33.*nx000.*nx011.*ny000+k22.*nx101.*ny000.^2+3.*k22.*nx000.*nx010.*ny001-k33.*nx000.*nx010.*ny001-k22.*nx100.*ny000.*ny001+k33.*nx100.*ny000.*ny001-k22.*nx000.*nx001.*ny010+k33.*nx000.*nx001.*ny010+2.*k33.*ny000.*ny001.*ny010-k11.*ny011+k22.*nx000.^2.*ny011+k33.*ny000.^2.*ny011+3.*k22.*nx001.*ny000.*ny100-k33.*nx001.*ny000.*ny100-2.*k22.*nx000.*ny001.*ny100+2.*k33.*nx000.*ny001.*ny100-k22.*nx000.*ny000.*ny101+k33.*nx000.*ny000.*ny101+k33.*nx001.^2.*nz000+2.*k22.*nx010.^2.*nz000-k33.*nx010.^2.*nz000+k22.*nx000.*nx020.*nz000-k33.*nx000.*nx020.*nz000-k22.*nx110.*ny000.*nz000+k33.*nx110.*ny000.*nz000+k33.*ny001.^2.*nz000-4.*k22.*nx010.*ny100.*nz000+2.*k33.*nx010.*ny100.*nz000+2.*k22.*ny100.^2.*nz000-k33.*ny100.^2.*nz000-k22.*nx000.*ny110.*nz000+k33.*nx000.*ny110.*nz000+k22.*ny000.*ny200.*nz000-k33.*ny000.*ny200.*nz000+k33.*nx101.*nz000.^2+k33.*ny011.*nz000.^2-k11.*nz002+(k24.*nz002)/2-2.*k22.*nx000.*nx010.*nz010+k22.*nx100.*ny000.*nz010-k33.*nx100.*ny000.*nz010-2.*k33.*ny000.*ny010.*nz010+k22.*nx000.*ny100.*nz010-k33.*nx000.*ny100.*nz010-k33.*nz000.*nz010.^2+(k24.*nz020)/2-k22.*nx000.^2.*nz020-k33.*ny000.^2.*nz020-k33.*nz000.^2.*nz020-2.*k33.*nx000.*nx100.*nz100+k22.*nx010.*ny000.*nz100-k33.*nx010.*ny000.*nz100+k22.*nx000.*ny010.*nz100-k33.*nx000.*ny010.*nz100-2.*k22.*ny000.*ny100.*nz100-k33.*nz000.*nz100.^2+2.*k22.*nx000.*ny000.*nz110-2.*k33.*nx000.*ny000.*nz110+(k24.*nz200)/2-k33.*nx000.^2.*nz200-k22.*ny000.^2.*nz200-k33.*nz000.^2.*nz200-2.*k22.*nx010.*q0+2.*k22.*ny100.*q0;

function [fex, fey, fez] = fFenn3D(nn,vv,dx,dy,dz,ea,e0)
%fFenn3D calculates the functional derivatives (used in frelax3D)
nx000 = nn(2:end-1,2:end-1,2:end-1,1);
ny000 = nn(2:end-1,2:end-1,2:end-1,2);
nz000 = nn(2:end-1,2:end-1,2:end-1,3);
v100 = (vv(4:end-1,3:end-2,3:end-2)-vv(2:end-3,3:end-2,3:end-2))/(2*dx);
v010 = (vv(3:end-2,4:end-1,3:end-2)-vv(3:end-2,2:end-3,3:end-2))/(2*dy);
v001 = (vv(3:end-2,3:end-2,4:end-1)-vv(3:end-2,3:end-2,2:end-3))/(2*dz);
% size(nz000),size(v001)

fex = e0.*ea.*nz000.*v001.*v100+e0.*ea.*ny000.*v010.*v100+e0.*ea.*nx000.*v100.^2;
fey = e0.*ea.*nz000.*v001.*v010+e0.*ea.*ny000.*v010.^2+e0.*ea.*nx000.*v010.*v100;
fez = e0.*ea.*nz000.*v001.^2+e0.*ea.*ny000.*v001.*v010+e0.*ea.*nx000.*v001.*v100;

function [vv] = fvoltageUpdate(nn,vv,dx,dy,dz,ea,eper)
%fvoltageUpdate updates voltage profile based on surrounding voltage and 
% director (used in frelax)
nx000 = nn(2:end-1,2:end-1,2:end-1,1);
ny000 = nn(2:end-1,2:end-1,2:end-1,2);
nz000 = nn(2:end-1,2:end-1,2:end-1,3);
nx100 = (nn(3:end,2:end-1,2:end-1,1)-nn(1:end-2,2:end-1,2:end-1,1))/(2*dx);
nx010 = (nn(2:end-1,3:end,2:end-1,1)-nn(2:end-1,1:end-2,2:end-1,1))/(2*dy);
nx001 = (nn(2:end-1,2:end-1,3:end,1)-nn(2:end-1,2:end-1,1:end-2,1))/(2*dz);
ny100 = (nn(3:end,2:end-1,2:end-1,2)-nn(1:end-2,2:end-1,2:end-1,2))/(2*dx);
ny010 = (nn(2:end-1,3:end,2:end-1,2)-nn(2:end-1,1:end-2,2:end-1,2))/(2*dy);
ny001 = (nn(2:end-1,2:end-1,3:end,2)-nn(2:end-1,2:end-1,1:end-2,2))/(2*dz);
nz100 = (nn(3:end,2:end-1,2:end-1,3)-nn(1:end-2,2:end-1,2:end-1,3))/(2*dx);
nz010 = (nn(2:end-1,3:end,2:end-1,3)-nn(2:end-1,1:end-2,2:end-1,3))/(2*dy);
nz001 = (nn(2:end-1,2:end-1,3:end,3)-nn(2:end-1,2:end-1,1:end-2,3))/(2*dz);

vr = vv(4:end-1,3:end-2,3:end-2);     vl = vv(2:end-3,3:end-2,3:end-2);
vf = vv(3:end-2,4:end-1,3:end-2);     vb = vv(3:end-2,2:end-3,3:end-2);
vu = vv(3:end-2,3:end-2,4:end-1);     vd = vv(3:end-2,3:end-2,2:end-3);

v100 = (vr-vl)/(2*dx);
v010 = (vf-vb)/(2*dy);
v001 = (vu-vd)/(2*dz);
v110 = (vv(4:end-1,4:end-1,3:end-2)-vv(4:end-1,2:end-3,3:end-2)-vv(2:end-3,4:end-1,3:end-2)+vv(2:end-3,2:end-3,3:end-2))/(4*dx*dy);
v101 = (vv(4:end-1,3:end-2,4:end-1)-vv(4:end-1,3:end-2,2:end-3)-vv(2:end-3,3:end-2,4:end-1)+vv(2:end-3,3:end-2,2:end-3))/(4*dx*dz);
v011 = (vv(3:end-2,4:end-1,4:end-1)-vv(3:end-2,4:end-1,2:end-3)-vv(3:end-2,2:end-3,4:end-1)+vv(3:end-2,2:end-3,2:end-3))/(4*dy*dz);

vv(3:end-2,3:end-2,3:end-2) = (eper*(dy^2*dz^2*(vl+vr)+dx^2*(dz^2*(vb+vf)+dy^2*(vd+vu)))+ea*(dy^2*dz^2*nx000.^2.*(vl+vr)+dx^2*(dz^2*ny000.^2.*(vb+vf)+dy^2*(dz^2*(nx100.*nz000.*v001+ny010.*nz000.*v001+2*nz000.*nz001.*v001+ny000.*nz010.*v001+nx000.*nz100.*v001+nx100.*ny000.*v010+2*ny000.*ny010.*v010+nx000.*ny100.*v010+ny001.*nz000.*v010+ny000.*nz001.*v010+2*ny000.*nz000.*v011+2*nx000.*nx100.*v100+nx010.*ny000.*v100+nx000.*ny010.*v100+nx001.*nz000.*v100+nx000.*nz001.*v100+2*nx000.*nz000.*v101+2*nx000.*ny000.*v110)+nz000.^2.*(vd+vu)))))./(2*(dy^2.*dz^2.*(eper+ea*nx000.^2)+dx^2.*(dz^2*(eper+ea*ny000.^2)+dy^2*(eper+ea*nz000.^2))));

function [fmx, fmy, fmz] = fFmnn3D(m,n,p,hx,hy,hz,u0)
%fFmnn3D calculates the functional derivatives (used in frelax3D)

fmx = ones(m-2,n-2,p-2)*-u0*hx*10000; % *10000 to convert to CGS: (G)
fmy = ones(m-2,n-2,p-2)*-u0*hy*10000;
fmz = ones(m-2,n-2,p-2)*-u0*hz*10000;

function [s, e, totalFE, fsk11, fsk22, fsk33, fsk24] = fFreeEnergy3D(hgui)
%fFreeEnergy3D calculates the elastic and electric free energy density
data = getappdata(hgui);

nn = data.nn;
vv = data.vv;

dx = str2double(data.edit_dimx)/str2double(data.edit_nx);
dy = str2double(data.edit_dimy)/str2double(data.edit_ny);
dz = str2double(data.edit_dimz)/str2double(data.edit_nz);

k11 = str2double(data.edit_k11);
k22 = str2double(data.edit_k22);
k33 = str2double(data.edit_k33);
k24 = str2double(data.edit_k24);

pitch = str2double(data.edit_pitch);
q0 = 2*pi/pitch;

epsper = str2double(data.edit_epsper);
ea = str2double(data.edit_epspar)-epsper;
e0 = 8.854e-12;

nx000 = nn(2:end-1,2:end-1,2:end-1,1);
ny000 = nn(2:end-1,2:end-1,2:end-1,2);
nz000 = nn(2:end-1,2:end-1,2:end-1,3);
nx100 = (nn(3:end,2:end-1,2:end-1,1)-nn(1:end-2,2:end-1,2:end-1,1))/(2*dx);
nx010 = (nn(2:end-1,3:end,2:end-1,1)-nn(2:end-1,1:end-2,2:end-1,1))/(2*dy);
nx001 = (nn(2:end-1,2:end-1,3:end,1)-nn(2:end-1,2:end-1,1:end-2,1))/(2*dz);
ny100 = (nn(3:end,2:end-1,2:end-1,2)-nn(1:end-2,2:end-1,2:end-1,2))/(2*dx);
ny010 = (nn(2:end-1,3:end,2:end-1,2)-nn(2:end-1,1:end-2,2:end-1,2))/(2*dy);
ny001 = (nn(2:end-1,2:end-1,3:end,2)-nn(2:end-1,2:end-1,1:end-2,2))/(2*dz);
nz100 = (nn(3:end,2:end-1,2:end-1,3)-nn(1:end-2,2:end-1,2:end-1,3))/(2*dx);
nz010 = (nn(2:end-1,3:end,2:end-1,3)-nn(2:end-1,1:end-2,2:end-1,3))/(2*dy);
nz001 = (nn(2:end-1,2:end-1,3:end,3)-nn(2:end-1,2:end-1,1:end-2,3))/(2*dz);
nx200 = (nn(3:end,2:end-1,2:end-1,1)+nn(1:end-2,2:end-1,2:end-1,1)-2*nx000)/dx^2;
nx020 = (nn(2:end-1,3:end,2:end-1,1)+nn(2:end-1,1:end-2,2:end-1,1)-2*nx000)/dy^2;
nx002 = (nn(2:end-1,2:end-1,3:end,1)+nn(2:end-1,2:end-1,1:end-2,1)-2*nx000)/dz^2;
ny200 = (nn(3:end,2:end-1,2:end-1,2)+nn(1:end-2,2:end-1,2:end-1,2)-2*ny000)/dx^2;
ny020 = (nn(2:end-1,3:end,2:end-1,2)+nn(2:end-1,1:end-2,2:end-1,2)-2*ny000)/dy^2;
ny002 = (nn(2:end-1,2:end-1,3:end,2)+nn(2:end-1,2:end-1,1:end-2,2)-2*ny000)/dz^2;
nz200 = (nn(3:end,2:end-1,2:end-1,3)+nn(1:end-2,2:end-1,2:end-1,3)-2*nz000)/dx^2;
nz020 = (nn(2:end-1,3:end,2:end-1,3)+nn(2:end-1,1:end-2,2:end-1,3)-2*nz000)/dy^2;
nz002 = (nn(2:end-1,2:end-1,3:end,3)+nn(2:end-1,2:end-1,1:end-2,3)-2*nz000)/dz^2;

fsk11 = (k11.*nx100.^2)/2+k11.*nx100.*ny010+(k11.*ny010.^2)/2+k11.*nx100.*nz001+k11.*ny010.*nz001+(k11.*nz001.^2)/2;
fsk22 = (k22.*nx001.^2.*ny000.^2)/2-k22.*nx000.*nx001.*ny000.*ny001+(k22.*nx000.^2.*ny001.^2)/2-k22.*nx001.*nx010.*ny000.*nz000+k22.*nx000.*nx010.*ny001.*nz000+k22.*nx001.*ny000.*ny100.*nz000-k22.*nx000.*ny001.*ny100.*nz000+(k22.*nx010.^2.*nz000.^2)/2-k22.*nx010.*ny100.*nz000.^2+(k22.*ny100.^2.*nz000.^2)/2+k22.*nx000.*nx001.*ny000.*nz010-k22.*nx000.^2.*ny001.*nz010-k22.*nx000.*nx010.*nz000.*nz010+k22.*nx000.*ny100.*nz000.*nz010+(k22.*nx000.^2.*nz010.^2)/2-k22.*nx001.*ny000.^2.*nz100+k22.*nx000.*ny000.*ny001.*nz100+k22.*nx010.*ny000.*nz000.*nz100-k22.*ny000.*ny100.*nz000.*nz100-k22.*nx000.*ny000.*nz010.*nz100+(k22.*ny000.^2.*nz100.^2)/2+k22.*nx001.*ny000.*q0-k22.*nx000.*ny001.*q0-k22.*nx010.*nz000.*q0+k22.*ny100.*nz000.*q0+k22.*nx000.*nz010.*q0-k22.*ny000.*nz100.*q0+(k22.*q0.^2)/2;
fsk33 = (k33.*nx000.^2.*nx001.^2)/2+(k33.*nx000.^2.*nx010.^2)/2+(k33.*nx010.^2.*ny000.^2)/2+k33.*nx000.*nx001.*ny000.*ny001+(k33.*ny000.^2.*ny001.^2)/2-k33.*nx000.^2.*nx010.*ny100-k33.*nx010.*ny000.^2.*ny100+(k33.*nx000.^2.*ny100.^2)/2+(k33.*ny000.^2.*ny100.^2)/2+k33.*nx001.*nx010.*ny000.*nz000-k33.*nx000.*nx010.*ny001.*nz000-k33.*nx001.*ny000.*ny100.*nz000+k33.*nx000.*ny001.*ny100.*nz000+(k33.*nx001.^2.*nz000.^2)/2+(k33.*ny001.^2.*nz000.^2)/2-k33.*nx000.*nx001.*ny000.*nz010-k33.*ny000.^2.*ny001.*nz010+k33.*nx000.*nx010.*nz000.*nz010-k33.*nx000.*ny100.*nz000.*nz010-k33.*ny001.*nz000.^2.*nz010+(k33.*ny000.^2.*nz010.^2)/2+(k33.*nz000.^2.*nz010.^2)/2-k33.*nx000.^2.*nx001.*nz100-k33.*nx000.*ny000.*ny001.*nz100-k33.*nx010.*ny000.*nz000.*nz100+k33.*ny000.*ny100.*nz000.*nz100-k33.*nx001.*nz000.^2.*nz100+k33.*nx000.*ny000.*nz010.*nz100+(k33.*nx000.^2.*nz100.^2)/2+(k33.*nz000.^2.*nz100.^2)/2;
fsk24 = (k24.*nx001.^2)/2+(k24.*nx000.*nx002)/2+(k24.*nx010.^2)/2+(k24.*nx000.*nx020)/2+(k24.*nx100.^2)/2+(k24.*nx000.*nx200)/2+(k24.*ny001.^2)/2+(k24.*ny000.*ny002)/2+k24.*nx100.*ny010+(k24.*ny010.^2)/2+(k24.*ny000.*ny020)/2-k24.*nx010.*ny100+(k24.*ny100.^2)/2+(k24.*ny000.*ny200)/2+k24.*nx100.*nz001+k24.*ny010.*nz001+(k24.*nz001.^2)/2+(k24.*nz000.*nz002)/2-k24.*ny001.*nz010+(k24.*nz010.^2)/2+(k24.*nz000.*nz020)/2-k24.*nx001.*nz100+(k24.*nz100.^2)/2+(k24.*nz000.*nz200)/2;

v100 = (vv(4:end-1,3:end-2,3:end-2)-vv(2:end-3,3:end-2,3:end-2))/(2*dx);
v010 = (vv(3:end-2,4:end-1,3:end-2)-vv(3:end-2,2:end-3,3:end-2))/(2*dy);
v001 = (vv(3:end-2,3:end-2,4:end-1)-vv(3:end-2,3:end-2,2:end-3))/(2*dz);

e = (e0.*epsper.*v001.^2)/2+(e0.*ea.*nz000.^2.*v001.^2)/2+e0.*ea.*ny000.*nz000.*v001.*v010+(e0.*epsper.*v010.^2)/2+(e0.*ea.*ny000.^2.*v010.^2)/2+e0.*ea.*nx000.*nz000.*v001.*v100+e0.*ea.*nx000.*ny000.*v010.*v100+(e0.*epsper.*v100.^2)/2+(e0.*ea.*nx000.^2.*v100.^2)/2;
s = fsk11+fsk22+fsk33-fsk24;

fs = -(k24.*nx001.^2)/2+(k33.*nx000.^2.*nx001.^2)/2-(k24.*nx000.*nx002)/2-(k24.*nx010.^2)/2+(k33.*nx000.^2.*nx010.^2)/2-(k24.*nx000.*nx020)/2+(k11.*nx100.^2)/2-(k24.*nx100.^2)/2-(k24.*nx000.*nx200)/2+(k22.*nx001.^2.*ny000.^2)/2+(k33.*nx010.^2.*ny000.^2)/2-k22.*nx000.*nx001.*ny000.*ny001+k33.*nx000.*nx001.*ny000.*ny001-(k24.*ny001.^2)/2+(k22.*nx000.^2.*ny001.^2)/2+(k33.*ny000.^2.*ny001.^2)/2-(k24.*ny000.*ny002)/2+k11.*nx100.*ny010-k24.*nx100.*ny010+(k11.*ny010.^2)/2-(k24.*ny010.^2)/2-(k24.*ny000.*ny020)/2+k24.*nx010.*ny100-k33.*nx000.^2.*nx010.*ny100-k33.*nx010.*ny000.^2.*ny100-(k24.*ny100.^2)/2+(k33.*nx000.^2.*ny100.^2)/2+(k33.*ny000.^2.*ny100.^2)/2-(k24.*ny000.*ny200)/2-k22.*nx001.*nx010.*ny000.*nz000+k33.*nx001.*nx010.*ny000.*nz000+k22.*nx000.*nx010.*ny001.*nz000-k33.*nx000.*nx010.*ny001.*nz000+k22.*nx001.*ny000.*ny100.*nz000-k33.*nx001.*ny000.*ny100.*nz000-k22.*nx000.*ny001.*ny100.*nz000+k33.*nx000.*ny001.*ny100.*nz000+(k33.*nx001.^2.*nz000.^2)/2+(k22.*nx010.^2.*nz000.^2)/2+(k33.*ny001.^2.*nz000.^2)/2-k22.*nx010.*ny100.*nz000.^2+(k22.*ny100.^2.*nz000.^2)/2+k11.*nx100.*nz001-k24.*nx100.*nz001+k11.*ny010.*nz001-k24.*ny010.*nz001+(k11.*nz001.^2)/2-(k24.*nz001.^2)/2-(k24.*nz000.*nz002)/2+k22.*nx000.*nx001.*ny000.*nz010-k33.*nx000.*nx001.*ny000.*nz010+k24.*ny001.*nz010-k22.*nx000.^2.*ny001.*nz010-k33.*ny000.^2.*ny001.*nz010-k22.*nx000.*nx010.*nz000.*nz010+k33.*nx000.*nx010.*nz000.*nz010+k22.*nx000.*ny100.*nz000.*nz010-k33.*nx000.*ny100.*nz000.*nz010-k33.*ny001.*nz000.^2.*nz010-(k24.*nz010.^2)/2+(k22.*nx000.^2.*nz010.^2)/2+(k33.*ny000.^2.*nz010.^2)/2+(k33.*nz000.^2.*nz010.^2)/2-(k24.*nz000.*nz020)/2+k24.*nx001.*nz100-k33.*nx000.^2.*nx001.*nz100-k22.*nx001.*ny000.^2.*nz100+k22.*nx000.*ny000.*ny001.*nz100-k33.*nx000.*ny000.*ny001.*nz100+k22.*nx010.*ny000.*nz000.*nz100-k33.*nx010.*ny000.*nz000.*nz100-k22.*ny000.*ny100.*nz000.*nz100+k33.*ny000.*ny100.*nz000.*nz100-k33.*nx001.*nz000.^2.*nz100-k22.*nx000.*ny000.*nz010.*nz100+k33.*nx000.*ny000.*nz010.*nz100-(k24.*nz100.^2)/2+(k33.*nx000.^2.*nz100.^2)/2+(k22.*ny000.^2.*nz100.^2)/2+(k33.*nz000.^2.*nz100.^2)/2-(k24.*nz000.*nz200)/2+k22.*nx001.*ny000.*q0-k22.*nx000.*ny001.*q0-k22.*nx010.*nz000.*q0+k22.*ny100.*nz000.*q0+k22.*nx000.*nz010.*q0-k22.*ny000.*nz100.*q0+(k22.*q0.^2)/2;
totalFE = fs-e;
% save('TFE','totalFE')

function c = pseudotensor(u,di)
% calculate the Handedness pseudotensor for a 3D director field

% Levi-Civita tensor
e = zeros(3,3,3);
e([8 12 22]) = 1;
e([6 16 20]) = -1;

shift = 0.5;
method = 'cubic';
[xdim,ydim,zdim,~] = size(u);
[Xq,Yq,Zq] = meshgrid(2:ydim-1,2:xdim-1,2:zdim-1);
upik(:,:,:,1,1) = interp3(squeeze(u(:,:,:,1)),Xq,Yq+shift,Zq,method);
umik(:,:,:,1,1) = interp3(squeeze(u(:,:,:,1)),Xq,Yq-shift,Zq,method);
upik(:,:,:,2,1) = interp3(squeeze(u(:,:,:,1)),Xq+shift,Yq,Zq,method);
umik(:,:,:,2,1) = interp3(squeeze(u(:,:,:,1)),Xq-shift,Yq,Zq,method);
upik(:,:,:,3,1) = interp3(squeeze(u(:,:,:,1)),Xq,Yq,Zq+shift,method);
umik(:,:,:,3,1) = interp3(squeeze(u(:,:,:,1)),Xq,Yq,Zq-shift,method);
upik(:,:,:,1,2) = interp3(squeeze(u(:,:,:,2)),Xq,Yq+shift,Zq,method);
umik(:,:,:,1,2) = interp3(squeeze(u(:,:,:,2)),Xq,Yq-shift,Zq,method);
upik(:,:,:,2,2) = interp3(squeeze(u(:,:,:,2)),Xq+shift,Yq,Zq,method);
umik(:,:,:,2,2) = interp3(squeeze(u(:,:,:,2)),Xq-shift,Yq,Zq,method);
upik(:,:,:,3,2) = interp3(squeeze(u(:,:,:,2)),Xq,Yq,Zq+shift,method);
umik(:,:,:,3,2) = interp3(squeeze(u(:,:,:,2)),Xq,Yq,Zq-shift,method);
upik(:,:,:,1,3) = interp3(squeeze(u(:,:,:,3)),Xq,Yq+shift,Zq,method);
umik(:,:,:,1,3) = interp3(squeeze(u(:,:,:,3)),Xq,Yq-shift,Zq,method);
upik(:,:,:,2,3) = interp3(squeeze(u(:,:,:,3)),Xq+shift,Yq,Zq,method);
umik(:,:,:,2,3) = interp3(squeeze(u(:,:,:,3)),Xq-shift,Yq,Zq,method);
upik(:,:,:,3,3) = interp3(squeeze(u(:,:,:,3)),Xq,Yq,Zq+shift,method);
umik(:,:,:,3,3) = interp3(squeeze(u(:,:,:,3)),Xq,Yq,Zq-shift,method);

diuk(:,:,:,1,1) = (upik(:,:,:,1,1)-umik(:,:,:,1,1))/2/di(1)/shift;
diuk(:,:,:,2,1) = (upik(:,:,:,2,1)-umik(:,:,:,2,1))/2/di(2)/shift;
diuk(:,:,:,3,1) = (upik(:,:,:,3,1)-umik(:,:,:,3,1))/2/di(3)/shift;
diuk(:,:,:,1,2) = (upik(:,:,:,1,2)-umik(:,:,:,1,2))/2/di(1)/shift;
diuk(:,:,:,2,2) = (upik(:,:,:,2,2)-umik(:,:,:,2,2))/2/di(2)/shift;
diuk(:,:,:,3,2) = (upik(:,:,:,3,2)-umik(:,:,:,3,2))/2/di(3)/shift;
diuk(:,:,:,1,3) = (upik(:,:,:,1,3)-umik(:,:,:,1,3))/2/di(1)/shift;
diuk(:,:,:,2,3) = (upik(:,:,:,2,3)-umik(:,:,:,2,3))/2/di(2)/shift;
diuk(:,:,:,3,3) = (upik(:,:,:,3,3)-umik(:,:,:,3,3))/2/di(3)/shift;

dimgrid = [xdim ydim zdim];
c = zeros([dimgrid 3 3]);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
c(2:end-1,2:end-1,2:end-1,i,j) = c(2:end-1,2:end-1,2:end-1,i,j)+...
    diuk(:,:,:,i,k)*e(j,l,k).*u(2:end-1,2:end-1,2:end-1,l);
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%% Callback functions for RELAX. %%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are called when the user interacts with the GUI

function edit_cd_Callback(hedit_cd,~,hgui)
% [edit_cd_Callback] Brief description:
    % Gets the string entered to the GUI. If cd string is empty reset the
    % cd string to the current folder using function: pwd. Try to reset the
    % current directory. If error set the current directory using pwd.
    % Update the edit_cd field the application data.

str = get(hedit_cd,'String');
if isempty(str) % if cd string is empty
    str = pwd;  % currentFolder = pwd; Identify the current folder
    set(hedit_cd, 'String', pwd);
end
try
    cd(str);
catch
    pwd;
    errordlg('Path not found','Warning');
    set(hedit_cd, 'String', pwd);
    return
end
setappdata(hgui,'edit_cd',pwd)    % update the application data
disp(['you changed directory to:  ' pwd])

function push_load_Callback(~,~,hgui)
% [push_load_Callback] Brief description:
    % 1.  Trys to load: 'filename' containing the variable: 'nn'. 
    % Requierments for nn: numel(size(nn)) = 4, m>=3||n>=3||p>=3||xyz123==3
    % 
    % 2.  Normalizes the director, calculates nx, ny, nz, dx, dy, dz and
    % the initial voltage profile, updates the application data, then runs
    % the plot Callback
    
% try to load the load filename with variable nn
try
	load(getappdata(hgui,'edit_load'),'nn');
catch
    % if this fails, display an error msg and return to the main function
	errordlg('File not found','Error');
	return
end

% if nn variable does not exist
if ~exist('nn','var')
    errordlg('nn not found','Error');
    return
end

% if nn variable is the wrong size
[m,n,p,xyz123]=size(nn); %#ok<NODEF>
if m<3 || n<3 || p<3 || xyz123~=3
    errordlg('nn wrong size, m<3 || n<3 || p<3 || xyz123~=3','Error');
    return
end

% Normalize the nn data from loaded file
mag = repmat(sqrt(sum((nn).^2,4)),1,1,1,3);
nn = nn./mag;
setappdata(hgui,'nn',nn)

% Update data.edit_n(i) with loaded array dimensions then calculate and
% update/create edit_dim(i) field: where (i) = {x,y,z}
% x
setappdata(hgui,'edit_nx',num2str(m))
setappdata(hgui,'edit_dx',...
    num2str(str2double(getappdata(hgui,'edit_dimx'))/m))
hx = findobj('tag','input_gridx');
set(hx,'String',num2str(m))
% y
setappdata(hgui,'edit_ny',num2str(n))
setappdata(hgui,'edit_dy',...
    num2str(str2double(getappdata(hgui,'edit_dimy'))/n))
hy = findobj('tag','input_gridy');
set(hy,'String',num2str(n))
% z
setappdata(hgui,'edit_nz',num2str(p))
setappdata(hgui,'edit_dz',...
    num2str(str2double(getappdata(hgui,'edit_dimz'))/p))
hz = findobj('tag','input_gridz');
set(hz,'String',num2str(p))

% Calculate the initial voltage profile fvoltageUpdate
u = str2double(getappdata(hgui,'edit_voltage'));
if u==0
    vv = zeros(m+2,n+2,p+2);
else
    dx = str2double(getappdata(hgui,'edit_dimx'))/m; % step size in um
    dy = str2double(getappdata(hgui,'edit_dimy'))/n;
    dz = str2double(getappdata(hgui,'edit_dimz'))/p;
    eper = str2double(getappdata(hgui,'edit_epsper'));%dielectric constants
    ea = str2double(getappdata(hgui,'edit_epspar'))-eper;
    [~,~,vv] = meshgrid(0:n+1,0:m+1,0:u/(p+1):u);
    vv = fvoltageUpdate(nn,vv,dx,dy,dz,ea,eper);
end
setappdata(hgui,'vv',vv)

% plot Callback to show director field after director is loaded
push_plot_Callback(1,1,hgui)

disp('you pushed load')
disp([num2str(m) 'X' num2str(n) 'X' num2str(p)])

function edit_load_Callback(hedit_load,~,hgui)
str = get(hedit_load,'String');
setappdata(hgui,'edit_load',str)
disp(['you entered string:  ' str])

function push_loadstate_Callback(~,~,hgui)

% try to load the state with data structure
try
	load(getappdata(hgui,'edit_loadstate'));
catch
    % if this fails, display an error msg and return to the main function
	errordlg('File not found','Error');
	return
end

% ----------------------------------------------------
%  tempararoly set children as inactive during process
children1 = findobj('tag','input','-or',...
    'tag','input_gridx','-or',...
    'tag','input_gridy','-or',...
    'tag','input_gridz','-or',...
    'tag','input_epspar','-or',...
    'tag','input_epsper','-or',...
    'tag','input_v','-or',...
    'tag','input_mx','-or',...
    'tag','input_my','-or',...
    'tag','input_mz','-or',...
    'tag','input_dp');
children2 = findobj('tag','input_val');
children3 = findobj('tag','input_k11','-or',...
    'tag','input_k22','-or',...
    'tag','input_k33','-or',...
    'tag','input_k24');
children4 = findobj('tag','input_dimx','-or',...
    'tag','input_dimy','-or',...
    'tag','input_dimz','-or',...
    'tag','input_pitch');
children = [children1;children2;children3;children4];
set(children,'Enable','inactive')
try
    for n = 1:length(children1)
        % get the callback name
        str = char(cell2mat(children1(n).Callback(1)));
        % clip callback
        str = ['data.' str(1:end-9)];
        % set GUI and appdata
        set(children1(n),'String',eval(str))
        setappdata(hgui,str(6:end),eval(str))
    end
    for n = 1:length(children2)
        % get the callback name
        str = char(cell2mat(children2(n).Callback(1)));
        % clip callback
        str = ['data.' str(1:end-9)];
        % set GUI and appdata
        set(children2(n),'Value',eval(str))
        setappdata(hgui,str(6:end),eval(str))
    end
    for n = 1:length(children3)
        % get the callback name
        str = char(cell2mat(children3(n).Callback(1)));
        % clip callback
        str = ['data.' str(1:end-9)];
        num = str2double(eval(str))*1e12; % 1e-12 -> 1
        % set GUI and appdata
        set(children3(n),'String',num2str(num))
        setappdata(hgui,str(6:end),num2str(num*1e-12)) % 1 -> 1e-12
    end
    for n = 1:length(children4)
        % get the callback name
        str = char(cell2mat(children4(n).Callback(1)));
        % clip callback
        str = ['data.' str(1:end-9)];
        num = str2double(eval(str))*1e6; % 1e-6 -> 1
        % set GUI and appdata
        set(children4(n),'String',num2str(num))
        setappdata(hgui,str(6:end),num2str(num*1e-6)) % 1 -> 1e-6
    end
catch
    errordlg('state file bad')
    set(children,'Enable','on')
    return
end

% hack to fix issue with load state----------------------------------------
% I fixed the cilderen3 loop and this hack may not be needed anymore,
% needes testing...
hpop_lc = findobj(children2(end));
pop_lc_Callback(hpop_lc,1,hgui)
hedit_k24 = findobj(hgui,'tag','input_k24');
edit_k24_Callback(hedit_k24,1,hgui)
% hack to fix issue with load state----------------------------------------

% set appdata (2)
if isfield(data,'edit_dx')
    nn = data.nn;
    vv = data.vv;
    setappdata(hgui,'nn',nn)
    setappdata(hgui,'vv',vv)
    dx = data.edit_dx;
    setappdata(hgui,'edit_dx',dx)
    dy = data.edit_dy;
    setappdata(hgui,'edit_dy',dy)
    dz = data.edit_dz;
    setappdata(hgui,'edit_dz',dz)
end

% set appearence (2)
hv = findobj('tag','input_v');
if str2double(hv(1).String)==0
    set(hv,'BackgroundColor',[1 1 1])
else
    set(hv,'BackgroundColor',[1 .7 .7])
end
hm = findobj('tag','input_mx','-or',...
    'tag','input_my','-or',...
    'tag','input_mz');
s1 = str2double(hm(3).String);
s2 = str2double(hm(2).String);
s3 = str2double(hm(1).String);
if s1==0
    set(hm(3),'BackgroundColor',[1 1 1])
else
    set(hm(3),'BackgroundColor',[.7 .7 1])
end
if s2==0
    set(hm(2),'BackgroundColor',[1 1 1])
else
    set(hm(2),'BackgroundColor',[.7 .7 1])
end
if s3==0
    set(hm(1),'BackgroundColor',[1 1 1])
else
    set(hm(1),'BackgroundColor',[.7 .7 1])
end

%  set children as on after process complete
set(children,'Enable','on')
% ----------------------------------------------------
if isfield(data,'nn') && ~isempty(data.nn)
    push_plot_Callback(1,1,hgui)
end

disp('you pushed load state')

function edit_loadstate_Callback(hedit_loadstate,~,hgui)
str = get(hedit_loadstate,'String');
setappdata(hgui,'edit_loadstate',str)
disp(['you entered string:  ' str])

function push_savestate_Callback(~,~,hgui)
data = getappdata(hgui); % make the data structure from all of the app data

if exist([data.edit_savestate '.mat'], 'file')
    % Construct a questdlg
    choice = questdlg('Would you like to overwrite?', ...
        'Options', ...
        'Yes','No','Cancel','Cancel');
    switch choice
        case 'Yes'
            save([data.edit_savestate '.mat'],'data')
        case 'No'
            warndlg('Must rename state (canceled)')
        case 'Cancel'
            return
    end
else
    save([data.edit_savestate '.mat'],'data')
end
disp('you pushed savestate')

function edit_savestate_Callback(hedit_savestate,~,hgui)
str = get(hedit_savestate,'String');
setappdata(hgui,'edit_savestate',str)
disp(['you entered string:  ' str])

function edit_nx_Callback(hedit_nx,~,hgui)
str = get(hedit_nx,'String');
setappdata(hgui,'edit_nx',str)
disp(['you entered string:  ' str])

function edit_ny_Callback(hedit_ny,~,hgui)
str = get(hedit_ny,'String');
setappdata(hgui,'edit_ny',str)
disp(['you entered string:  ' str])

function edit_nz_Callback(hedit_nz,~,hgui)
str = get(hedit_nz,'String');
setappdata(hgui,'edit_nz',str)
disp(['you entered string:  ' str])

function edit_dimx_Callback(hedit_dimx,~,hgui)
str = get(hedit_dimx,'String');
setappdata(hgui,'edit_dimx',num2str(str2double(str)*1e-6))
setappdata(hgui,'edit_dx',num2str(str2double(str)*1e-6/...
    str2double(getappdata(hgui,'edit_nx'))))
disp(['you entered string:  ' str])

function edit_dimy_Callback(hedit_dimy,~,hgui)
str = get(hedit_dimy,'String');
setappdata(hgui,'edit_dimy',num2str(str2double(str)*1e-6))
setappdata(hgui,'edit_dy',num2str(str2double(str)*1e-6/...
    str2double(getappdata(hgui,'edit_ny'))))
disp(['you entered string:  ' str])

function edit_dimz_Callback(hedit_dimz,~,hgui)
d = get(hedit_dimz,'String');
setappdata(hgui,'edit_dimz',num2str(str2double(d)*1e-6))
setappdata(hgui,'edit_dz',num2str(str2double(d)*1e-6/...
    str2double(getappdata(hgui,'edit_nz'))))
h1 = findobj('tag','input_pitch');
p = str2double(get(h1,'String'));
dp = str2double(d)/p;
setappdata(hgui,'edit_dp',num2str(dp))
h2 = findobj('tag','input_dp');
set(h2,'String',num2str(dp))
disp(['you entered string:  ' d])

function pop_lc_Callback(hpop_lc,~,hgui)
val = get(hpop_lc,'Value');
children = findobj('tag','input_k11','-or',...
                   'tag','input_k22','-or',...
                   'tag','input_k33','-or',...
                   'tag','input_epspar','-or',...
                   'tag','input_epsper');
switch val
    case 1 % 5CB
        k11 = 6.4e-12; % splay elastic constant [pN]
        k22 = 3.0e-12; % twist Frank elastic constant [pN]
        k33 = 10.e-12; % bend Frank elastic constant [pN]
        set(children(5),'String',num2str(k11*10e11))
        set(children(4),'String',num2str(k22*10e11))
        set(children(3),'String',num2str(k33*10e11))
        setappdata(hgui,'edit_k11',num2str(k11))
        setappdata(hgui,'edit_k22',num2str(k22))
        setappdata(hgui,'edit_k33',num2str(k33))
        epar = 19.0; % relative dielectric permittivity para[dimensionless]
        eper = 5.2; % relative dielectric permittivity perp [dimensionless]
        set(children(1),'String',num2str(eper))
        set(children(2),'String',num2str(epar))
        setappdata(hgui,'edit_epsper',num2str(eper))
        setappdata(hgui,'edit_epspar',num2str(epar))
    case 2 % ZLI2806
        k11 = 14.9e-12; % splay elastic constant [pN]
        k22 = 7.9e-12; % twist Frank elastic constant [pN]
        k33 = 15.4e-12; % bend Frank elastic constant [pN]
        set(children(5),'String',num2str(k11*10e11))
        set(children(4),'String',num2str(k22*10e11))
        set(children(3),'String',num2str(k33*10e11))
        setappdata(hgui,'edit_k11',num2str(k11))
        setappdata(hgui,'edit_k22',num2str(k22))
        setappdata(hgui,'edit_k33',num2str(k33))
        epar = 3.3; % relative dielectric permittivity para [dimensionless]
        eper = 8.1; % relative dielectric permittivity perp [dimensionless]
        set(children(1),'String',num2str(eper))
        set(children(2),'String',num2str(epar))
        setappdata(hgui,'edit_epsper',num2str(eper))
        setappdata(hgui,'edit_epspar',num2str(epar))
    case 3 % AMLC0010
        k11 = 17.2e-12; % splay elastic constant [pN]
        k22 = 7.51e-12; % twist Frank elastic constant [pN]
        k33 = 17.9e-12; % bend Frank elastic constant [pN]
        set(children(5),'String',num2str(k11*10e11))
        set(children(4),'String',num2str(k22*10e11))
        set(children(3),'String',num2str(k33*10e11))
        setappdata(hgui,'edit_k11',num2str(k11))
        setappdata(hgui,'edit_k22',num2str(k22))
        setappdata(hgui,'edit_k33',num2str(k33))
        epar = 3.4; % relative dielectric permittivity para [dimensionless]
        eper = 7.1; % relative dielectric permittivity perp [dimensionless]
        set(children(1),'String',num2str(eper))
        set(children(2),'String',num2str(epar))
        setappdata(hgui,'edit_epsper',num2str(eper))
        setappdata(hgui,'edit_epspar',num2str(epar))
    case 4 % User Defined
        k11 = str2double(get(children(5),'String'))*10e-13; 
        k22 = str2double(get(children(4),'String'))*10e-13; 
        k33 = str2double(get(children(3),'String'))*10e-13; 
        setappdata(hgui,'edit_k11',num2str(k11))
        setappdata(hgui,'edit_k22',num2str(k22))
        setappdata(hgui,'edit_k33',num2str(k33))
end
setappdata(hgui,'pop_lc',val)
disp(['you entered LC value:  ' num2str(val)])

function edit_k11_Callback(hedit_k11,~,hgui)
str = num2str(str2double(get(hedit_k11,'String'))*10e-13);
setappdata(hgui,'edit_k11',str)
disp(['you entered string:  ' str])

function edit_k22_Callback(hedit_k22,~,hgui)
str = num2str(str2double(get(hedit_k22,'String'))*10e-13);
setappdata(hgui,'edit_k22',str)
disp(['you entered string:  ' str])

function edit_k33_Callback(hedit_k33,~,hgui)
str = num2str(str2double(get(hedit_k33,'String'))*10e-13);
setappdata(hgui,'edit_k33',str)
disp(['you entered string:  ' str])

function edit_k24_Callback(hedit_k24,~,hgui)
str = num2str(str2double(get(hedit_k24,'String'))*10e-13);
setappdata(hgui,'edit_k24',str)
disp(['you entered string:  ' str])

function edit_pitch_Callback(hedit_pitch,~,hgui)
p = str2double(get(hedit_pitch,'String'))*1e-6;
setappdata(hgui,'edit_pitch',num2str(p))
h1 = findobj('tag','input_dimz');
d = str2double(get(h1,'String'))*1e-6;
dp = d/p;
setappdata(hgui,'edit_dp',num2str(dp))
h2 = findobj('tag','input_dp');
set(h2,'String',num2str(dp))
disp(['you entered string:  ' num2str(p*1e6)])

function edit_dp_Callback(hedit_dp,~,hgui)
dp = str2double(get(hedit_dp,'String'));
setappdata(hgui,'edit_dp',num2str(dp))
h1 = findobj('tag','input_dimz');
d = str2double(get(h1,'String'));
p = d/dp;
setappdata(hgui,'edit_pitch',num2str(p*1e-6))
h2 = findobj('tag','input_pitch');
set(h2,'String',num2str(p))
disp(['you entered string:  ' num2str(dp)])

function edit_epspar_Callback(hedit_epspar,~,hgui)
str = get(hedit_epspar,'String');
setappdata(hgui,'edit_epspar',str)
disp(['you entered string:  ' str])

function edit_epsper_Callback(hedit_epsper,~,hgui)
str = get(hedit_epsper,'String');
setappdata(hgui,'edit_epsper',str)
disp(['you entered string:  ' str])

function edit_voltage_Callback(hedit_voltage,~,hgui)

str = get(hedit_voltage,'String');
setappdata(hgui,'edit_voltage',str)
u = str2double(str);
if u~=0
    set(hedit_voltage,'BackgroundColor',[1 .7 .7])
else
    set(hedit_voltage,'BackgroundColor',[1 1 1])
end

% Calculate the initial voltage profile using fvoltageUpdate
nn = getappdata(hgui,'nn');
[m,n,p,~]=size(nn);
if u==0
    vv = zeros(m+2,n+2,p+2);
else
    dx = str2double(getappdata(hgui,'edit_dimx'))/m; % step size in um
    dy = str2double(getappdata(hgui,'edit_dimy'))/n;
    dz = str2double(getappdata(hgui,'edit_dimz'))/p;
    eper = str2double(getappdata(hgui,'edit_epsper'));%dielectric constants
    ea = str2double(getappdata(hgui,'edit_epspar'))-eper;
    [~,~,vv] = meshgrid(0:n+1,0:m+1,0:u/(p+1):u);
    vv = fvoltageUpdate(nn,vv,dx,dy,dz,ea,eper);
end
setappdata(hgui,'vv',vv)

disp(['you entered string:  ' str])

function edit_mx_Callback(hedit_mx,~,hgui)
str = get(hedit_mx,'String');
setappdata(hgui,'edit_mx',str)
if str2double(str)~=0
    set(hedit_mx,'BackgroundColor',[.7 .7 1])
else
    set(hedit_mx,'BackgroundColor',[1 1 1])
end
disp(['you entered string:  ' str])

function edit_my_Callback(hedit_my,~,hgui)
str = get(hedit_my,'String');
setappdata(hgui,'edit_my',str)
if str2double(str)~=0
    set(hedit_my,'BackgroundColor',[.7 .7 1])
else
    set(hedit_my,'BackgroundColor',[1 1 1])
end
disp(['you entered string:  ' str])

function edit_mz_Callback(hedit_mz,~,hgui)
str = get(hedit_mz,'String');
setappdata(hgui,'edit_mz',str)
if str2double(str)~=0
    set(hedit_mz,'BackgroundColor',[.7 .7 1])
else
    set(hedit_mz,'BackgroundColor',[1 1 1])
end
disp(['you entered string:  ' str])

function rb_periodicBCsTB_Callback(hrb_periodicBCsTB,~,hgui)
test = (get(hrb_periodicBCsTB,'Value')== get(hrb_periodicBCsTB,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_periodicBCsTB',test)

function rb_periodicBCsLR_Callback(hrb_periodicBCsLR,~,hgui)
test = (get(hrb_periodicBCsLR,'Value')== get(hrb_periodicBCsLR,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_periodicBCsLR',test)

function rb_periodicBCsFB_Callback(hrb_periodicBCsFB,~,hgui)
test = (get(hrb_periodicBCsFB,'Value')== get(hrb_periodicBCsFB,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_periodicBCsFB',test)

function push_plot_Callback(~,~,hgui)
% [push_plot_Callback] Brief description:
    % Makes the status text BackgroundColor red while plotting then invokes
    % the plotter and after turns tatus text BackgroundColor back to blue

htext_status = findobj('tag','htext_status');
set(htext_status,'BackgroundColor','red','String','Busy')
drawnow % ensures that the color apears during plotting
ploter(hgui) % function to plot director
set(htext_status,'BackgroundColor','blue','String','Ready')

function pop_style_Callback(hpop_style,~,hgui)
val = get(hpop_style,'Value');
setappdata(hgui,'pop_style',val)
disp(['you entered style value:  ' num2str(val)])

function pop_xsection_Callback(hpop_xsection,~,hgui)
val = get(hpop_xsection,'Value');
setappdata(hgui,'pop_xsection',val)
disp(['you entered xsection value:  ' num2str(val)])

function edit_slicenum_Callback(hedit_slicenum,~,hgui)
str = get(hedit_slicenum,'String');
setappdata(hgui,'edit_slicenum',str)
disp(['you entered string:  ' str])

function edit_down_sample_Callback(hedit_down_sample,~,hgui)
str = get(hedit_down_sample,'String');
set(hedit_down_sample,'String',str)
setappdata(hgui,'edit_down_sample',str)
disp(['you entered string:  ' str])

function pop_energy_type_Callback(hpop_energy_type,~,hgui)
val = get(hpop_energy_type,'Value');
setappdata(hgui,'pop_energy_type',val)
disp(['you entered xsection value:  ' num2str(val)])

function edit_interp_energy_Callback(hedit_interp_energy,~,hgui)
str = get(hedit_interp_energy,'String');
setappdata(hgui,'edit_interp_energy',str)
disp(['you entered string:  ' str])

function rb_cutoff_Callback(hrb_cutoff,~,hgui)
test = (get(hrb_cutoff,'Value')== get(hrb_cutoff,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_cutoff',test)

function edit_cutoff_Callback(hedit_cutoff,~,hgui)
str = get(hedit_cutoff,'String');
if str2double(str)<=0
    str = getappdata(hgui,'edit_cutoff');
    set(hedit_cutoff,'String',str)
    errordlg('cutoff must be > 0')
end
setappdata(hgui,'edit_cutoff',str)
disp(['you entered string:  ' str])

function push_print_Callback(~,~,hgui)
str = getappdata(hgui,'edit_print');
style = getappdata(hgui,'pop_style');
FEon = getappdata(hgui,'pop_energy_type');

if style==5 || FEon>1 % 3D vectors or plot FE
    renderer = '-opengl';
else
    renderer = '-painters';
end
print('-depsc','-tiff','-r300',renderer,str)
disp(['you pushed print:  ' str])

function edit_print_Callback(hedit_print,~,hgui)
str = get(hedit_print,'String');
setappdata(hgui,'edit_print',str)
disp(['you entered string:  ' str])

function edit_iterations_Callback(hedit_iterations,~,hgui)
str = get(hedit_iterations,'String');
setappdata(hgui,'edit_iterations',str)
disp(['you entered string:  ' str])

function rb_ms_Callback(hrb_ms,~,hgui)
test = (get(hrb_ms,'Value')== get(hrb_ms,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_ms',test)

function rb_txtmsg_Callback(hrb_txtmsg,~,hgui)
test = (get(hrb_txtmsg,'Value')== get(hrb_txtmsg,'Max'));
if  test
    display('Selected');
else
    display('Not selected');
end
setappdata(hgui,'rb_txtmsg',test)

function edit_gain_Callback(hedit_gain,~,hgui)
str = get(hedit_gain,'String');
setappdata(hgui,'edit_gain',str)
disp(['you entered string:  ' str])

function push_relax_Callback(~,~,hgui)
tic
htext_it = findobj('tag','text_it');
set(htext_it,'String',1)
htext_total_it = findobj('tag','htext_total_it');
set(htext_total_it,'String',getappdata(hgui,'edit_iterations'))
htext_status = findobj('tag','htext_status');
set(htext_status,'BackgroundColor','red','String','Busy')
drawnow

% relax function here
num_it = str2double(getappdata(hgui,'edit_iterations'));
num_rep =  1;
flag = 1;
kk = 0;
while flag==1
    for report = 1:num_rep
        % relax the data
        frelax(hgui,report,num_it,num_rep);
        % update reports
        str = num_it/num_rep + (report-1)*num_it/num_rep;
        set(htext_it,'String',str)
        % plot the data
        ploter(hgui)
    end % report = 1:num_rep
    flag = getappdata(hgui,'rb_txtmsg');
% %     save([getappdata(hgui,'data.edit_savestate') '.mat'],'data')
    kk = kk+1;
    disp(['Continuous mode: ' num2str(kk)])
end % while flag==1

set(htext_status,'BackgroundColor','blue','String','Ready')
% elapsed time output
elapsedTime = toc;
if elapsedTime > 3600
    elapsedTime = [num2str(toc/3600) '(h)'];
elseif  elapsedTime > 60
    elapsedTime = [num2str(toc/60) '(m)'];
else
    elapsedTime = [num2str(toc) '(s)'];
end % elapsedTime
htext_elapsed_time = findobj('tag','text_elapsed_time');
set(htext_elapsed_time,'String',['Elapsed time: ' elapsedTime])

% test internet connection
[~,flag]=urlread('http://www.google.com');
% if radio button txtmsg, and connected to internet, send text to number
if getappdata(hgui,'rb_txtmsg')==1 && flag==1
    txtmsg('6092735335','verizon','Relax',['Number of iterations:  ' ... 
       'status: DONE'])
elseif flag~=1 && getappdata(hgui,'rb_txtmsg')==1
    errordlg('Internet Connection Error','Error');
end % data.rb_txtmsg==1 && flag==1
disp('you pushed RELAX')

%% %%%%%%%%%%%%%%%%%%% Plotting functions for RELAX. %%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are used for plotting

function ploter(hgui)
data = getappdata(hgui);

[m,n,p,dim] = size(data.nn);
if m<3 || n<3 || p<3 || dim <3
    errordlg('must load director (nn at least 3X3X3X3)')
    return
end

% sample and slice the data
[dataaspect,colorx,colory,maxnum,w,v,u,x,y] = slicer(hgui);

% set up the axis
ha = getappdata(hgui,'axis1');
axes(ha),cla,colorbar off
set(ha,'xcolor',colorx,'ycolor',colory,...
    'dataaspectratio',dataaspect);

hold on
% plot energy or handedness if within the data
if (1<str2double(data.edit_slicenum))&&...
        (str2double(data.edit_slicenum)<maxnum)
    plotIMAGE(hgui)
end
% plot the director field
plotdirector(hgui,w,v,u,x,y)
hold off

% fix axes position and limits
M = str2double(data.edit_down_sample);
switch data.pop_xsection
    case 1
        if m>n
            axis([-1,max([m*M n*M]),...
                -1+diff([m*M n*M])/2,max([m*M n*M])+diff([m*M n*M])/2])
        else
            axis([-1-diff([m*M n*M])/2,max([m*M n*M])-diff([m*M n*M])/2,...
                -1,max([m*M n*M])])
        end
    case 2
        if m>p
            axis([-1,max([m*M p*M]),...
                -1+diff([m*M p*M])/2,max([m*M p*M])+diff([m*M p*M])/2])
        else
            axis([-1-diff([m*M p*M])/2,max([m*M p*M])-diff([m*M p*M])/2,...
                -1,max([m*M p*M])])
        end
    case 3
        if n>p
            axis([-1,max([n*M p*M]),...
                -1+diff([n*M p*M])/2,max([n*M p*M])+diff([n*M p*M])/2])
        else
            axis([-1-diff([n*M p*M])/2,max([n*M p*M])-diff([n*M p*M])/2,...
                -1,max([n*M p*M])])
        end
end
zoom reset % this makes doubleclick zoom function work for new plot limits
drawnow

function [dataaspect,colorx,colory,maxnum,w,v,u,x,y] = slicer(hgui)
data = getappdata(hgui);
% get plot parameters
xsection = data.pop_xsection;
slicenum = str2double(data.edit_slicenum);
rate = str2double(data.edit_down_sample);
aspectx = str2double(data.edit_dimx)/str2double(data.edit_nx);
aspecty = str2double(data.edit_dimy)/str2double(data.edit_ny);
aspectz = str2double(data.edit_dimz)/str2double(data.edit_nz);
% sample nn by rate
nn(:,:,:,1) = data.nn(:,:,:,1);
nn(:,:,:,2) = data.nn(:,:,:,2);
nn(:,:,:,3) = data.nn(:,:,:,3);
[m,n,p,~] = size(nn);

% create data cross section
switch xsection
    case 1
        dataaspect = [aspecty aspectx 1];
        xx = linspace(1,m,m*rate);
        yy = linspace(1,n,n*rate);
        [xx,yy,zz] = meshgrid(yy,xx,slicenum);
        % make the sampled nn slice
        w = squeeze(interp3(nn(:,:,:,3),xx,yy,zz,'cubic'));
        v = squeeze(interp3(nn(:,:,:,1),xx,yy,zz,'cubic'));
        u = squeeze(interp3(nn(:,:,:,2),xx,yy,zz,'cubic'));
        [xdim,ydim] = size(v);
        [x,y] = meshgrid(1:ydim,1:xdim);
        colorx = [0 0 1];
        colory = [1 0 0];
        maxnum = p;
    case 2
        dataaspect = [aspectz aspectx 1];
        xx = linspace(1,m,m*rate);
        zz = linspace(1,p,p*rate);
        [xx,yy,zz] = meshgrid(slicenum,xx,zz);
        % make the sampled nn slice
        w = -squeeze(interp3(nn(:,:,:,2),xx,yy,zz,'cubic'));
        v = squeeze(interp3(nn(:,:,:,1),xx,yy,zz,'cubic'));
        u = squeeze(interp3(nn(:,:,:,3),xx,yy,zz,'cubic'));
        [xdim,ydim] = size(v);
        [x,y] = meshgrid(1:ydim,1:xdim);
        colorx = [0 0 1];
        colory = [0 1 0];
        maxnum = n;
    case 3
        dataaspect = [aspecty aspectz 1];
        yy = linspace(1,n,n*rate);
        zz = linspace(1,p,p*rate);
        [xx,yy,zz] = meshgrid(yy,slicenum,zz);
        % make the sampled nn slice
        w = squeeze(interp3(nn(:,:,:,1),xx,yy,zz,'cubic'));
        v = squeeze(interp3(nn(:,:,:,2),xx,yy,zz,'cubic'));
        u = squeeze(interp3(nn(:,:,:,3),xx,yy,zz,'cubic'));
        [xdim,ydim] = size(v);
        [x,y] = meshgrid(1:ydim,1:xdim);
        colorx = [1 0 0];
        colory = [0 1 0];
        maxnum = m;
end
x = x-1;
y = y-1;

function plotIMAGE(hgui)

plot_type = getappdata(hgui,'pop_energy_type'); % type of plot
switch plot_type
    case 1
        return
    case 2
        [~, ~, totalFE, ~, ~, ~, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,totalFE)
    case 3
        [~, ~, ~, fsk11, ~, ~, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,fsk11)
    case 4
        [~, ~, ~, ~, fsk22, ~, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,fsk22)
    case 5
        [~, ~, ~, ~, ~, fsk33, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,fsk33)
    case 6
        [~, ~, ~, ~, ~, ~, fsk24] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,-fsk24)
    case 7
        [s, ~, ~, ~, ~, ~, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,s)
    case 8
        [~, e, ~, ~, ~, ~, ~] = fFreeEnergy3D(hgui);
        % function to plot energy with colorbar
        plotEnergy(hgui,-e)
    case 9
        dx = str2double(getappdata(hgui,'edit_dimx'))/...
            str2double(getappdata(hgui,'edit_nx'));
        dy = str2double(getappdata(hgui,'edit_dimy'))/...
            str2double(getappdata(hgui,'edit_ny'));
        dz = str2double(getappdata(hgui,'edit_dimz'))/...
            str2double(getappdata(hgui,'edit_nz'));
        di = [dx dy dz];
        c = pseudotensor(getappdata(hgui,'nn'),di);
        % function to plot Handedness pseudotensor with colorbar cool
        plotHandedness(hgui,c)
    otherwise
        return
end

function plotHandedness(hgui,c)
data = getappdata(hgui); % load the data for ".notation"

slicenum = str2double(data.edit_slicenum);
M = str2double(getappdata(hgui,'edit_down_sample'));
% % M = str2double(data.edit_down_sample);

c11 = c(:,:,:,1,1);
c22 = c(:,:,:,2,2);
c33 = c(:,:,:,3,3);
tr = c11+c22+c33;
[m,n,p] = size(tr);

% make the sampled handedness pseudotensor trace slice
switch data.pop_xsection
    case 1
        xx = linspace(2,m-1,(m-2));
        yy = linspace(2,n-1,(n-2));
        [xx,yy,zz] = meshgrid(yy,xx,slicenum);
        trslice = squeeze(interp3(tr,xx,yy,zz,'cubic'))';
    case 2
        xx = linspace(2,m-1,(m-2));
        zz = linspace(2,p-1,(p-2));
        [xx,yy,zz] = meshgrid(slicenum,xx,zz);
        trslice = squeeze(interp3(tr,xx,yy,zz,'cubic'))';
    case 3
        yy = linspace(2,n-1,(n-2));
        zz = linspace(2,p-1,(p-2));
        [xx,yy,zz] = meshgrid(yy,slicenum,zz);
        trslice = squeeze(interp3(tr,xx,yy,zz,'cubic'))';
end

[y,x] = size(trslice);
interp = str2double(data.edit_interp_energy);
if x==1 || y==1
    [Xq,Yq] = meshgrid(1:x,1:y);
else
    [Xq,Yq] = meshgrid(1:1/interp:x,1:1/interp:y);
    trslice = interp2(trslice,Xq,Yq,'cubic');
end

% % if max(trslice(:))<2e6 || min(trslice(:))>-2e6
% %     cutoff = min(abs([max(trslice(:)) min(trslice(:))]));
% %     trslice(abs(trslice)>=cutoff) = cutoff;
% % else
% %     cutoffhigh = 2e6;
% %     cutofflow = -2e6;
% %     trslice(trslice>=cutoffhigh)=cutoffhigh;
% %     trslice(trslice<=cutofflow)=cutofflow;
% % end
% % caxis('auto')

pitch = str2double(data.edit_pitch);
caxis([-3*2*pi/abs(pitch) 3*2*pi/abs(pitch)])
map = flipud(cool);
colormap(map)
imagesc(Xq(1,:)*M+(M-1)*0.5,Yq(:,1)*M+(M-1)*0.5,trslice)
hcb=colorbar;
set(hcb,'YTick',-3*2*pi/abs(pitch):2*pi/abs(pitch):3*2*pi/abs(pitch))
set(hcb,'YTickLabel',{'-3|p|' '-2|p|' '-|p|' '0' '|p|' '2|p|' '3|p|'})

function plotEnergy(hgui,energy)

data = getappdata(hgui); % load the data for ".notation"

slicenum = str2double(data.edit_slicenum);
[m,n,p] = size(energy);
switch data.pop_xsection
    case  1 % 'XY'
        xx = linspace(1,m,m);
        yy = linspace(1,n,n);
        [yy,xx,zz] = meshgrid(yy,xx,slicenum);
        [x,y,z] = size(energy);
        if x==1 || y==1
            return
        elseif z==1
            energySlice = squeeze(interp2(squeeze(energy),squeeze(yy)',squeeze(xx)','cubic'));
        else
            energySlice = squeeze(interp3(energy,yy,xx,zz,'cubic'))';
        end
    case  2 % 'XZ'
        xx = linspace(1,m,m);
        zz = linspace(1,p,p);
        [yy,xx,zz] = meshgrid(slicenum,xx,zz);
        [x,y,z] = size(energy);
        if x==1 || z==1
            return
        elseif y==1
            energySlice = squeeze(interp2(squeeze(energy),squeeze(zz)',squeeze(xx)','cubic'));
        else
            energySlice = squeeze(interp3(energy,yy,xx,zz,'cubic'))';
        end
    case  3 % 'YZ'
        xx = linspace(1,n,n);
        zz = linspace(1,p,p);
        [yy,xx,zz] = meshgrid(slicenum,xx,zz);
        [x,y,z] = size(energy);
        if y==1 || z==1
            return
        elseif x==1
            energySlice = squeeze(interp2(squeeze(energy),squeeze(zz)',squeeze(yy)','cubic'));
        else
            energySlice = squeeze(interp3(energy,yy,xx,zz,'cubic'))';
        end
end

interp = str2double(data.edit_interp_energy);
[y,x]=size(energySlice);
Xq = linspace(1,x,x*interp);
Yq = linspace(1,y,y*interp);
[Xq,Yq] = meshgrid(Xq,Yq);
energySlice = interp2(energySlice,Xq,Yq,'cubic');
% adjust the parula colormap to make variations in energy near 0 more
% visable and avoid red or blue so cylinders caps are still visable
map = parula+...
    repmat([0.0 0.3 0.0],length(parula),1)-...
    repmat([0.1 0.0 0.3],length(parula),1); 
rb_cutoff = data.rb_cutoff; % saturate color at this value 
if rb_cutoff==1
    cutoff = str2double(data.edit_cutoff); % saturate color at this value
    energySlice(abs(energySlice)>cutoff) = cutoff;
    if max(energySlice(:))==cutoff
        map(end,:)=[1 1 1];
    end
end
map(map<0) = 0;
map(map>1) = 1;
caxis([min(energySlice(:))-2*eps max(energySlice(:))])%-2*eps prevents err
colormap(map)
M = str2double(getappdata(hgui,'edit_down_sample'));
imagesc(Xq(1,:)*M+(M-1)*0.5,Yq(:,1)*M+(M-1)*0.5,energySlice)
colorbar

function plotdirector(hgui,w,v,u,x,y)
data = getappdata(hgui);
% get plot parameters
style = data.pop_style;
% plot director
switch style
    case 1 % stick
        quiver(y-v/2,x-u/2,v,u,'k.','autoscale','off');
    case 2 % stick/dot
        plot(y,x,'k.')
        quiver(y-v/2,x-u/2,v,u,'k.','autoscale','off');
    case 3 % nail
        plot(y(w>0)+v(w>0)/2,x(w>0)+u(w>0)/2,'ro',...
            'MarkerSize',4,'MarkerFaceColor',[1 0.5 0.5])
        plot(y(w<0)-v(w<0)/2,x(w<0)-u(w<0)/2,'bo',...
            'MarkerSize',4,'MarkerFaceColor',[0.5 0.5 1])
        quiver(y-v/2,x-u/2,v,u,'k.','autoscale','off');

    case 4 % vector
        plot(y-v/2,x-u/2,'k.')
        quiver(y-v/2,x-u/2,v,u,'k','autoscale','off');

    case 5 % 3Dvector
        DDDvectors(w,v,u)
        
    case 6 % cylinder
        cylinders(w,v,u)
        
    case 7
        disp('none')
end

function DDDvectors(w,v,u)
U = v;
V = u;
W = w;
arrow_color = zeros(numel(W),3);
for ind=1:numel(V)
    if V(ind)>=0
        arrow_color(ind,:) = [1 0 0];
    elseif V(ind)<0
        arrow_color(ind,:) = [0 0 1];
    end
end
scaling = 1.4;
N = 12;
dim = size(U');
[Y,X,Z] = meshgrid(1:dim(1),1:dim(2),1);
quiver3d(X-1-U/2,Y-1-V/2,Z-W/2,U,V+eps,W,arrow_color,scaling,N);  % +eps for some cases!!!!!!!!!!!
light('Style','local','Position',[1 0 1],'Color',[1 0 1]);
material([.7 1 .7]) 

function cylinders(w,v,u)
% [cylinders] Brief description:
%   Input director data to plot cylinders with red and blue ends using
%   MATLABs patch. This function uses vectorization for speed.

r = 0.12; % radius of cylinder
res = 8; % at least 1.5
c3 = [.7 .7 .7]; % color of body (grey)
theta = -pi:pi/res:pi;

nn(:,:,1) = v';
nn(:,:,2) = w';
nn(:,:,3) = u';

[x,y,~] = size(nn);

% make color data for cylinder caps
red = nn(:,:,2)>0;   % Red if oriented out of the page
blue = nn(:,:,2)<=0; % Blue if oriented into the page
color = zeros(2*x*y,3); % initialize cdata true color -> [red green blue]
color(1:end/2,1) = red(:); % red if oriented out of the page 
color(end/2+1:end,3) = blue(:); % blue if oriented into the page
color(color(:,1)==0&color(:,2)==0&color(:,3)==0,:)=.7; % otherwize grey
% sort the caps so grey prints first and does not cover the colored caps
[~,reorder] = sort(sum(color,2),'descend'); 
color = color(reorder,:); % this reorers the cdata

p1 = cat(3, nn(:,:,1)/2, nn(:,:,2)/2, nn(:,:,3)/2);    % center of cap1
p2 = cat(3,-nn(:,:,1)/2,-nn(:,:,2)/2,-nn(:,:,3)/2);    % center of cap2

% make a normalized vector basis
p(:,:,3) = reshape(abs(nn(:,:,1))>0.7,size(nn(:,:,1)));
p(:,:,1) = reshape(abs(nn(:,:,2))>0.7,size(nn(:,:,2)));
p(:,:,2) = sum(p,3)==0; % a point used to define the vector basis
R1 = cross(p-p1,p2-p1,3)./...
    reshape(repmat(sqrt(sum(cross(p-p1,p2-p1,3).^2,3)),1,3),size(nn));
S1 = cross(R1,p2-p1,3)./...
    reshape(repmat(sqrt(sum(cross(R1,p2-p1,3).^2,3)),1,3),size(nn));

% define the cylinder body (rectangles)
p1x = p1(:,:,1);
p1y = p1(:,:,3);
p2x = p2(:,:,1);
p2y = p2(:,:,3);

[xx,yy] = meshgrid(1:y,1:x);

% create 3 rectangles to make the cylinder body
a1 = cos(fatan(p1x,p1y)+pi/2)*r;
b1 = sin(fatan(p1x,p1y)+pi/2)*r;
a2 = cos(fatan(p2x,p2y)+pi/2)*r;
b2 = sin(fatan(p2x,p2y)+pi/2)*r;
rect1x = shiftdim(reshape([xx+p1x+a1,xx+p1x+a1,xx+p2x-a2,xx+p2x+a2]-1,...
    [x,y,4]),2);
rect1y = shiftdim(reshape([yy+p1y+b1,yy+p1y+b1,yy+p2y-b2,yy+p2y+b2]-1,...
    [x,y,4]),2);
rect2x = shiftdim(reshape([xx+p1x+a1,xx+p1x-a1,xx+p2x+a2,xx+p2x-a2]-1,...
    [x,y,4]),2);
rect2y = shiftdim(reshape([yy+p1y+b1,yy+p1y-b1,yy+p2y+b2,yy+p2y-b2]-1,...
    [x,y,4]),2);
rect3x = shiftdim(reshape([xx+p1x-a1,xx+p1x-a1,xx+p2x+a2,xx+p2x-a2]-1,...
    [x,y,4]),2);
rect3y = shiftdim(reshape([yy+p1y-b1,yy+p1y-b1,yy+p2y+b2,yy+p2y-b2]-1,...
    [x,y,4]),2);
rectx = [rect1x(:);rect2x(:);rect3x(:)];
recty = [rect1y(:);rect2y(:);rect3y(:)];

% combine the vertices and faces for fast patch
verts =  [rectx(:) recty(:)];
faces = reshape(1:length(verts),4,x*y*3)';
patchinfo.Vertices = verts;
patchinfo.Faces = faces;
patchinfo.FaceColor = c3;
patch(patchinfo,'EdgeColor','None');

% create a basis that is the right size for calculation of cap vertices
r1x = reshape(repmat(R1(:,:,1),1,1,length(theta)),[x y length(theta)]);
r1y = reshape(repmat(R1(:,:,3),1,1,length(theta)),[x y length(theta)]);
s1x = reshape(repmat(S1(:,:,1),1,1,length(theta)),[x y length(theta)]);
s1y = reshape(repmat(S1(:,:,3),1,1,length(theta)),[x y length(theta)]);
th = reshape(repmat(theta,x*y,1),[x,y,length(theta)]);
[xx,yy] = meshgrid(1:y,1:x,1:length(theta));

% define the red cap vertices
p1x = reshape(repmat(p1(:,:,1),1,1,length(theta)),[x y length(theta)]);
p1y = reshape(repmat(p1(:,:,3),1,1,length(theta)),[x y length(theta)]);
Qx1 = shiftdim(xx+p1x+r*cos(th).*r1x+r*sin(th).*s1x-1,2);
Qy1 = shiftdim(yy+p1y+r*cos(th).*r1y+r*sin(th).*s1y-1,2);

% define the blue cap vertices
p2x = reshape(repmat(p2(:,:,1),1,1,length(theta)),[x y length(theta)]);
p2y = reshape(repmat(p2(:,:,3),1,1,length(theta)),[x y length(theta)]);
Qx2 = shiftdim(xx+p2x+r*cos(th).*r1x+r*sin(th).*s1x-1,2);
Qy2 = shiftdim(yy+p2y+r*cos(th).*r1y+r*sin(th).*s1y-1,2);

% combine the vertices and faces of all caps
vx = [Qx1(:);Qx2(:)];
vy = [Qy1(:);Qy2(:)];
verts = [vx(:) vy(:)];
faces = reshape(1:length(verts),length(theta),2*x*y)';
patchcaps.Vertices = verts;
% reorder to print grey caps then colored caps
patchcaps.Faces = faces(reorder,:);
patch(patchcaps,...
    'FaceVertexCData',color,...
    'FaceColor','flat',...
    'EdgeColor','None'); 

%% %%%%%%%%%%%%%%%%%% Subroutine functions for RELAX. %%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are used by other functions

function A=fatan(X,Y)
% fatan.m
% proper arctangent!
% daa 081020
a=(X<0);  % 0 in I, IV; 1 in II, III
b=(-1).^a;  % 1 in I, IV; -1 in II, III
c=(-1).^(Y<0);  % 1 in I, II; -1 in III, IV
A=-pi*a.*b.*c+atan(Y./X);

function txtmsg(number,carrier,subject,message)
% txtmsg:  example >> txtmsg('6092735335','verizon','test','message')
% Note:  All inputs are string format. Preferences are for GMAIL accounts.
%   Input:  number (10 digit string w/o dashes spaces etc.)
%           carrier (select your carrier) info from ->
% http://www.sms411.net/2006/07/how-to-send-email-to-phone.html
%           subject (optional) 
%           message 
% txtmsg('6092735335','verizon','Relax',['txt' 'status: DONE'])
% =========================================================================
% YOUR EMAIL AND PASSWORD HERE: (your account connected to your phone)
myaddress = 'smalyukhlabnotifications@gmail.com';
mypassword = 'tictorons';
% =========================================================================
% test internet connection
[~,flag]=urlread('http://www.google.com');
% if connected to internet send text to number
if flag~=1
    errordlg('Internet Connection Error','Error');
end

if nargin == 3 % if no subject (only 3 inputs)
    message = subject;
    subject = '';
end

switch strrep(strrep(lower(carrier),'-',''),'&','')
    case 'alltel';    emailto = strcat(number,'@message.alltel.com');
    case 'att';       emailto = strcat(number,'@mmode.com');
    case 'boost';     emailto = strcat(number,'@myboostmobile.com');
    case 'cingular';  emailto = strcat(number,'@cingularme.com');
    case 'cingular2'; emailto = strcat(number,'@mobile.mycingular.com');
    case 'nextel';    emailto = strcat(number,'@messaging.nextel.com');
    case 'sprint';    emailto = strcat(number,'@messaging.sprintpcs.com');
    case 'tmobile';   emailto = strcat(number,'@tmomail.net');
    case 'verizon';   emailto = strcat(number,'@vtext.com');
    case 'virgin';    emailto = strcat(number,'@vmobl.com');
end

% Set up the preferences:
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);
% For using gmail
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send the text mesage
sendmail(emailto,subject,message)

function p=quiver3d(X,Y,Z,U,V,W,arrow_color,scaling,N)

% Linearize the input data...
  X=X(:); Y=Y(:); Z=Z(:);
  U=U(:); V=V(:); W=W(:);

% Input color can be single color vector (1x3 or 3x1) *OR* indexed color
% vector (1xN or Nx1, where N=length(X)) *OR* true color matrix (3xN or Nx3)
% Regardless, the input color has to be repeated for each arrow

single_color =isequal([3 1],size(arrow_color(:)));                                                                                         % Single  color therefore [3x1] or [1x3]...
indexed_color=(size(arrow_color,1)==1 & size(arrow_color,2)==length(X(:))) | (size(arrow_color,1)==length(X(:)) & size(arrow_color,2)==1); % Indexed color therefore [1xN] or [Nx1]...
true_color   =(size(arrow_color,1)==3 & size(arrow_color,2)==length(X(:))) | (size(arrow_color,1)==length(X(:)) & size(arrow_color,2)==3); % True    color therefore [3xN] or [Nx3]...

if single_color
  arrow_color=repmat(arrow_color(:)',[size(X,1) 1]);
elseif indexed_color
  arrow_color=arrow_color(:);                                 
elseif true_color 
  if     isequal([3 3],size(arrow_color)), warning('Ambiguous TRUE COLOR definition intended results *may* require transpose of color matrix.');
  elseif size(arrow_color,1)==3, arrow_color=arrow_color';
  end
else 
  error('Color argument did not appear to be a SINGLE color, nor INDEXED color, nor TRUE color.');
end

[xat,yat,zat,faces0,~]=gen_template_arrow(X,Y,Z,N);

D=0.5*mean((U.^2+V.^2+W.^2).^0.5);   % Calculate an arrow body diameter base on the mean vector magnitude

vertices=[]; faces=[]; fc=[];

for k=1:size(X,1)

    % Scale the normalized arrow data by the vector's length...then autoscale the data
    xa=xat*sqrt(U(k,1).^2+V(k,1).^2+W(k,1).^2);
    ya=yat*D;
    za=zat*D;

    if scaling
      A=get_autoscale_value(X,Y,Z,U,V,W,scaling);
      xa=A*xa;
      ya=A*ya;
      za=A*za;
    end      

  % Generate orthogonal basis for the rotation
    Evct(:,1)=[U(k,1); V(k,1); W(k,1)]/norm([U(k,1); V(k,1); W(k,1)]);  % First unit vector along vector direction
    Evct(:,2)=cross(Evct(:,1),[1 0 0])/norm(cross(Evct(:,1),[1 0 0]));
    Evct(:,3)=cross(Evct(:,1),Evct(:,2));

  % Rotate the template arrow
    XYZ=Evct*[xa(:)'; ya(:)'; za(:)'];
    [xa,ya,za]=deal(reshape(XYZ(1,:),size(xa)),reshape(XYZ(2,:),size(ya)),reshape(XYZ(3,:),size(za)));
  
  % Translate the template arrow
    xa=xa+X(k);  ya=ya+Y(k);  za=za+Z(k);  % x,y,z are the tesselated surface points...

  % Triangulate the surface points
  vertices0=[xa(:),ya(:),za(:)];   

  fc0=repmat(arrow_color(k,:),[size(faces0,1) 1]);  
  
  % Concatenate the patch surfaces for each glyph
  vertices=[vertices; vertices0                     ]; %#ok<AGROW>
  faces   =[faces;    faces0+(k-1)*size(vertices0,1)]; %#ok<AGROW>
  fc      =[fc;       fc0                           ]; %#ok<AGROW>
  
end

% Draw all the glyphs
  p=patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',fc,'FaceColor','Flat');
    set(p,'EdgeColor','None');
    set(p,'BackFaceLighting','lit','FaceLighting','Phong');
    set(p,'SpecularColorReflectance',0.1,'SpecularExponent',1);
    set(p,'DiffuseStrength',0.75);

function [xa,ya,za,faces,vertices] = gen_template_arrow(~,~,~,N)
% Generate a normalized arrow geometry
% L=1;             % Arrow length...
% R=1/10;          % Arrow radius
% S=3/10;          % Arrow slope
% 
% L=3/4;             % Arrow length...
% R=3/40;            % Arrow radius...
% S=1/4;             % Arrow slope...

L=.26;             % Arrow length...
R=.07;             % Arrow radius...
S=.25;             % Arrow slope...

% Generate the FAUX surface data for tesselation.
% ...can't easily tesselate the true surface data because of the
% ...sudden change in radius at the Tip-Body interface
[xt,yt,zt]=cylinder(linspace(R,0,2),max(N));  % Arrow Tip
[xb,yb,zb]=cylinder([R 2*R],max(N));          % Arrow Body
zb=zb*L;                                      % Scale the body
zt=S*zt+2*L;                                  % Scale and displace the arrow Tip

xx=[zt zb];                                    % Combine the body and top coordinates
yy=[yt yb];                                    % Combine the body and top coordinates
zz=[xt xb];                                    % Combine the body and top coordinates

vertices=[xx(:),yy(:),zz(:)];
faces=convhulln(vertices);

% Generate the REAL surface data for rendering
[xt,yt,zt]=cylinder(linspace(0.7*S,0,2),max(N));  % Arrow tip
[xb,yb,zb]=cylinder(R,max(N));                    % Arrow body

% Scale the data and shift the tip to the end of the arrow body...
zb=zb*L;
zt=S*zt+L;

xa=[zt zb];  % Final arrow coordinates (tip+body)
ya=[yt yb];
za=[xt xb];

% The get_auto_scale function uses code that was borrowed from The
% Mathwork's QUIVER3 function.
function autoscale = get_autoscale_value(x,y,z,u,v,w,autoscale)
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n; 
  dely = diff([min(y(:)) max(y(:))])/m;
  delz = diff([min(z(:)) max(y(:))])/max(m,n);
  del = sqrt(delx.^2 + dely.^2 + delz.^2);
  if del>0
    len = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
    maxlen = max(len(:));
  else
    maxlen = 0;
  end
  
  if maxlen>0
    autoscale = autoscale*0.9 / maxlen;
  else
    autoscale = autoscale*0.9;
  end

%% %%%%%%%%%%%%%%%%%%%% Generate the GUI for RELAX. %%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions generate the GUI when the relax function is first run

function hgui = make_fig_fcn
% make_fig_fcn Brief description:
%   Generates the GUI axis and the GUI organization. Initializes the data
%   structure in the aplication data. Outputs the GUI handle

% figure and figure properties
hgui = figure('Visible','off'); % generate an invisible figure
hgui.IntegerHandle = 'off'; % turn off integer handle to prevent overwrite
hgui.Color = [0.87 0.9 0.93]; % background color
hgui.DockControls = 'off'; % prevents GUI from "docking" into matlab
hgui.MenuBar = 'figure'; % turn off the default menubar
hgui.Name = mfilename; % set the GUI name to file name "relax"
hgui.Units = 'normalized';
hgui.OuterPosition = [0 0 1 1]; % [left bottom width height]
hgui.Units = 'pixels';
gui.OuterPosition.full = hgui.OuterPosition; % gui size (full screen)
gui.OuterPosition.default = [1 1 940 800]; % gui size (default)
sizeflag = gui.OuterPosition.default>gui.OuterPosition.full; % test display
if sum(sizeflag)==0 % if defalt size fits on screen
    hgui.OuterPosition = gui.OuterPosition.default;
else % if defalt size is too big for screen
    for lbwh = 1:4
        if sizeflag(lbwh)~=0
            gui.OuterPosition.default(lbwh) = gui.OuterPosition.full(lbwh);
        end
    end
end
gui.lbwh = hgui.OuterPosition; % measured 'pixels' gui position
hgui.Resize = 'on';

% initialize data structure
timenow = datetime('now','TimeZone','local','Format','yMMdd_HHmmss');
disp(mfilename)
disp(timenow)
setappdata(hgui,'Name',mfilename)     % Function file name
setappdata(hgui,'Time',char(timenow)) % DateTime (YYYYMMDD_hhmmss)
setappdata(hgui,'nn',[])              % initialize nn
setappdata(hgui,'vv',[])              % initialize vv
setappdata(hgui,'gui',gui)            % gui information

% Axis #1
w = 750; % axis dimensions (in pixles)
left = 175/gui.lbwh(3);  % normalized left position
bottom = 15/gui.lbwh(4); % normalized bottom position
width = w/gui.lbwh(3);   % normalized width
height = w/gui.lbwh(4);  % normalized height
axis1 = axes('Parent',hgui,...
    'Position',[left bottom width height],...
    'tag','axis',...
    'color','none',...
    'box','on','xcolor',[0 0 1],'ycolor',[1 0 0],...
    'linewidth',1,...
    'ticklength',[0 0],'xticklabel','','yticklabel','');
setappdata(hgui,'axis1',axis1) % axis1 information

% Axis #2 (plot ms parameter)
axis2 = axes('Parent',hgui,...
    'OuterPosition',[.01 .01 .18 .18],...
    'tag','axis_ms');
setappdata(hgui,'axis2',axis2) % axis2 information

elmt = 0;
% current directory
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = cd;
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [1 1 20 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_cd_Callback,hgui};
setappdata(hgui,'edit_cd',ui(elmt).str)
% current directory lable
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'Current Directory';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [1.2 21.1 4 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% load pushbutton
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'load';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [2 1 2 1]; % row column norm_width norm_hight
ui(elmt).callback = {@push_load_Callback,hgui};
% load filename input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = 'filename';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [2 3.1 4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_load_Callback,hgui};
setappdata(hgui,'edit_load',ui(elmt).str)

% loadstate pushbutton
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'load state';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [3 1 2 1]; % row column norm_width norm_hight
ui(elmt).callback = {@push_loadstate_Callback,hgui};
% state filename input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = 'state';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [3 3.1 4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_loadstate_Callback,hgui};
setappdata(hgui,'edit_loadstate',ui(elmt).str)

% savestate pushbutton
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'save state';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [4 1 2 1]; % row column norm_width norm_hight
ui(elmt).callback = {@push_savestate_Callback,hgui};
% savestate filename input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = 'state';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [4 3.1 4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_savestate_Callback,hgui};
setappdata(hgui,'edit_savestate',ui(elmt).str)

% grid label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'size of computational grid';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [5.2 1 6 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% gridx num input
elmt = elmt+1;
ui(elmt).tag = 'input_gridx';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [6 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_nx_Callback,hgui};
setappdata(hgui,'edit_nx',ui(elmt).str)
% gridx num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'nx';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [6.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% gridy num input
elmt = elmt+1;
ui(elmt).tag = 'input_gridy';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [6 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_ny_Callback,hgui};
setappdata(hgui,'edit_ny',ui(elmt).str)
% gridy num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'ny';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [6.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% gridz num input
elmt = elmt+1;
ui(elmt).tag = 'input_gridz';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [6 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_nz_Callback,hgui};
setappdata(hgui,'edit_nz',ui(elmt).str)
% gridz num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'nz';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [6.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% dimx input
elmt = elmt+1;
ui(elmt).tag = 'input_dimx';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [7 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_dimx_Callback,hgui};
setappdata(hgui,'edit_dimx','1e-6')
% gridx num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'um';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [7.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% gridy num input
elmt = elmt+1;
ui(elmt).tag = 'input_dimy';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [7 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_dimy_Callback,hgui};
setappdata(hgui,'edit_dimy','1e-6')
% gridy num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'um';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [7.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% gridz num input
elmt = elmt+1;
ui(elmt).tag = 'input_dimz';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [7 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_dimz_Callback,hgui};
setappdata(hgui,'edit_dimz','1e-6')
% gridz num label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'um';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [7.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% LC pop
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'popupmenu';
ui(elmt).str = {'5CB','ZLI2806','AMLC0010','UserDefined'};
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [8 1.3 5 1]; % row column norm_width norm_hight
ui(elmt).callback = {@pop_lc_Callback,hgui};
ui(elmt).value = 4;
setappdata(hgui,'pop_lc',ui(elmt).value)

% k11 input
elmt = elmt+1;
ui(elmt).tag = 'input_k11';
ui(elmt).style = 'edit';
ui(elmt).str = '10';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [9 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_k11_Callback,hgui};
setappdata(hgui,'edit_k11',num2str(str2double(ui(elmt).str)*10e-13))
% k11 label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'k11';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [9.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% k22 input
elmt = elmt+1;
ui(elmt).tag = 'input_k22';
ui(elmt).style = 'edit';
ui(elmt).str = '10';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [9 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_k22_Callback,hgui};
setappdata(hgui,'edit_k22',num2str(str2double(ui(elmt).str)*10e-13))
% k22 label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'k22';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [9.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% k33 input
elmt = elmt+1;
ui(elmt).tag = 'input_k33';
ui(elmt).style = 'edit';
ui(elmt).str = '10';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [9 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_k33_Callback,hgui};
setappdata(hgui,'edit_k33',num2str(str2double(ui(elmt).str)*10e-13))
% k33 label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'k33';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [9.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% k24 input
elmt = elmt+1;
ui(elmt).tag = 'input_k24';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [10 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_k24_Callback,hgui};
setappdata(hgui,'edit_k24',num2str(str2double(ui(elmt).str)*10e-13))
% k24 label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'k24';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [10.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% pitch input
elmt = elmt+1;
ui(elmt).tag = 'input_pitch';
ui(elmt).style = 'edit';
ui(elmt).str = 'inf';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [10 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_pitch_Callback,hgui};
setappdata(hgui,'edit_pitch',ui(elmt).str)
% pitch label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'p';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [10.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% dp input
elmt = elmt+1;
ui(elmt).tag = 'input_dp';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [10 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_dp_Callback,hgui};
setappdata(hgui,'edit_dp',ui(elmt).str)
% dp label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'd/p';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [10.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';
  
% epspar input
elmt = elmt+1;
ui(elmt).tag = 'input_epspar';
ui(elmt).style = 'edit';
ui(elmt).str = '10';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [11 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_epspar_Callback,hgui};
setappdata(hgui,'edit_epspar',ui(elmt).str)
% epspar label
elmt = elmt+1;
ui(elmt).tag = 'symbol';
ui(elmt).style = 'text';
ui(elmt).str = 'e||';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [11.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% epsper input
elmt = elmt+1;
ui(elmt).tag = 'input_epsper';
ui(elmt).style = 'edit';
ui(elmt).str = '5';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [11 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_epsper_Callback,hgui};
setappdata(hgui,'edit_epsper',ui(elmt).str)
% epsper label
elmt = elmt+1;
ui(elmt).tag = 'symbol';
ui(elmt).style = 'text';
ui(elmt).str = 'e^';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [11.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% voltage input
elmt = elmt+1;
ui(elmt).tag = 'input_v';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [11 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_voltage_Callback,hgui};
setappdata(hgui,'edit_voltage',ui(elmt).str)
% voltage label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'V';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [11.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% mx input
elmt = elmt+1;
ui(elmt).tag = 'input_mx';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [12 1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_mx_Callback,hgui};
setappdata(hgui,'edit_mx',ui(elmt).str)
% mx label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'mx';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [12.1 2 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% my input
elmt = elmt+1;
ui(elmt).tag = 'input_my';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [12 3 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_my_Callback,hgui};
setappdata(hgui,'edit_my',ui(elmt).str)
% my label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'my';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [12.1 4 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% mz input
elmt = elmt+1;
ui(elmt).tag = 'input_mz';
ui(elmt).style = 'edit';
ui(elmt).str = '0';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [12 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_mz_Callback,hgui};
setappdata(hgui,'edit_mz',ui(elmt).str)
% mz label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'mz';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [12.1 6 1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% periodicBCs label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'i/o periodic BCs';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [13.2 1 6 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% periodicBCsTB radio
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'TB';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [14 1 1.5 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_periodicBCsTB_Callback,hgui};
ui(elmt).value = 0;
setappdata(hgui,'rb_periodicBCsTB',ui(elmt).value)

% periodicBCsLR radio
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'LR';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [14 3 1.5 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_periodicBCsLR_Callback,hgui};
ui(elmt).value = 0;
setappdata(hgui,'rb_periodicBCsLR',ui(elmt).value)

% periodicBCsFB radio
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'FB';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [14 5 1.5 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_periodicBCsFB_Callback,hgui};
ui(elmt).value = 0;
setappdata(hgui,'rb_periodicBCsFB',ui(elmt).value)

% plot pushbutton
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'plot';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [15 1 2 1]; % row column norm_width norm_hight
ui(elmt).callback = {@push_plot_Callback,hgui};

% style pop
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'popupmenu';
ui(elmt).str = {'stick','stick/dot','nail','vector','3Dvector',...
    'cylinder','none'};
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [15 3.1 4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@pop_style_Callback,hgui};
ui(elmt).value = 1;
setappdata(hgui,'pop_style',ui(elmt).value)

% xsection pop
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'popupmenu';
ui(elmt).str = {'XY','XZ','YZ'};
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [16 .8 2.485 1]; % row column norm_width norm_hight
ui(elmt).callback = {@pop_xsection_Callback,hgui};
ui(elmt).value = 1;
setappdata(hgui,'pop_xsection',ui(elmt).value)

% slicenum input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [16 3.1 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_slicenum_Callback,hgui};
setappdata(hgui,'edit_slicenum',ui(elmt).str)
% slicenum label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = '#';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [16.1 4.15 .9 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% down_sample input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [16 5 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_down_sample_Callback,hgui};
setappdata(hgui,'edit_down_sample',ui(elmt).str)
% down_sample label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'M';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [16.1 6 .9 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% pop_energy_type pop
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'popupmenu';
ui(elmt).str = {'none','total FE','k11','k22','k33','k24','elastic',...
    'electric','handedness'};
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [17 .8 3.4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@pop_energy_type_Callback,hgui};
ui(elmt).value = 1;
setappdata(hgui,'pop_energy_type',ui(elmt).value)

% edit_interp_energy input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [17 4.05 1 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_interp_energy_Callback,hgui};
setappdata(hgui,'edit_interp_energy',ui(elmt).str)
% interp_energy label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'interp';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [17.1 5.15 2 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% rb_cutoff radio
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'cutoff';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [18 1 3.4 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_cutoff_Callback,hgui};
ui(elmt).value = 0;
setappdata(hgui,'rb_cutoff',ui(elmt).value)

% edit_cutoff input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '10';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [18 3.5 3.5 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_cutoff_Callback,hgui};
setappdata(hgui,'edit_cutoff',ui(elmt).str)

% print pushbutton
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'print';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [19 1 2 1]; % row column norm_width norm_hight
ui(elmt).callback = {@push_print_Callback,hgui};

% print input
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = 'printname';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [19 3.1 3.9 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_print_Callback,hgui};
setappdata(hgui,'edit_print',ui(elmt).str)

% iterations label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'iterations';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [21.3 1 3 1]; % row column norm_width norm_hight
ui(elmt).callback = '';
% edit_iterations
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [22 1 2.9 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_iterations_Callback,hgui};
setappdata(hgui,'edit_iterations',ui(elmt).str)

% edit_gain label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'MSTS Gain';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [21.3 4 3.1 1]; % row column norm_width norm_hight
ui(elmt).callback = '';
% edit_gain
elmt = elmt+1;
ui(elmt).tag = 'input';
ui(elmt).style = 'edit';
ui(elmt).str = '1';
ui(elmt).bgc = [.94 .94 .94];
ui(elmt).rcwh = [22 4 3 1]; % row column norm_width norm_hight
ui(elmt).callback = {@edit_gain_Callback,hgui};
setappdata(hgui,'edit_gain',ui(elmt).str)

% push_relax pushbutton
elmt = elmt+1;
ui(elmt).tag = 'relax';
ui(elmt).style = 'pushbutton';
ui(elmt).str = 'RELAX';
ui(elmt).bgc = 'green';
ui(elmt).rcwh = [24 4 3 2]; % row column norm_width norm_hight
ui(elmt).callback = {@push_relax_Callback,hgui};

% htext_it label
elmt = elmt+1;
ui(elmt).tag = 'text_it';
ui(elmt).style = 'text';
ui(elmt).str = '1';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [23 1 3 1]; % row column norm_width norm_hight
ui(elmt).callback = '';
% text_of label
elmt = elmt+1;
ui(elmt).tag = '';
ui(elmt).style = 'text';
ui(elmt).str = 'of';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [23.6 1 3 1]; % row column norm_width norm_hight
ui(elmt).callback = '';
% htext_total_it label
elmt = elmt+1;
ui(elmt).tag = 'htext_total_it';
ui(elmt).style = 'text';
ui(elmt).str = '1';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [24.2 1 3 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

% text_status
elmt = elmt+1;
ui(elmt).tag = 'htext_status';
ui(elmt).style = 'text';
ui(elmt).str = 'Ready';
ui(elmt).bgc = 'blue';
ui(elmt).rcwh = [26 1 3 2]; % row column norm_width norm_hight
ui(elmt).callback = {@text_status_Callback,hgui};

% rb_ms radio
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'Ms';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [25 4 3 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_ms_Callback,hgui};
ui(elmt).value = 1;
setappdata(hgui,'rb_ms',ui(elmt).value)

% rb_txtmsg
elmt = elmt+1;
ui(elmt).tag = 'input_val';
ui(elmt).style = 'radiobutton';
ui(elmt).str = 'continuous';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [26 4 3 1]; % row column norm_width norm_hight
ui(elmt).callback = {@rb_txtmsg_Callback,hgui};
ui(elmt).value = 0;
setappdata(hgui,'rb_txtmsg',ui(elmt).value)

% text_elapsed_time
elmt = elmt+1;
ui(elmt).tag = 'text_elapsed_time';
ui(elmt).style = 'text';
ui(elmt).str = 'Elapsed time: (s)';
ui(elmt).bgc = hgui.Color;
ui(elmt).rcwh = [27 1 6 1]; % row column norm_width norm_hight
ui(elmt).callback = '';

make_uicontrols(hgui,ui)

function make_uicontrols(hgui,ui)
% [make_uicontrols] Brief description:
%   Generates the GUI UIcontrols from make_fig_fcn organization
% 
gui = getappdata(hgui,'gui');
leftmar = 4;
rowhight = 18;
columnwidth = 25;
for n=1:length(ui)
    topmar = 104+2*(ui(n).rcwh(1)-1);
    row = ui(n).rcwh(1);
    column = ui(n).rcwh(2);
    columns = ui(n).rcwh(3);
    rows = ui(n).rcwh(4);
    left = (column-1)*columnwidth+leftmar;
    bottom = gui.lbwh(4)-(row-1)*rowhight-topmar-rowhight;
    width = columnwidth*columns;
    height = rowhight*rows;
    pos = [left bottom width height];
    if isempty(ui(n).value) && ~strcmp(ui(n).tag,'symbol')
        uicontrol(hgui,...
        'tag',ui(n).tag,...
        'Style',ui(n).style,...
        'String',ui(n).str,...
        'BackgroundColor',ui(n).bgc,...
        'Position',pos,...
        'Callback',ui(n).callback);
    elseif strcmp(ui(n).tag,'symbol')
        uicontrol(hgui,...
        'tag',ui(n).tag,...
        'Style',ui(n).style,...
        'String',ui(n).str,'FontName','symbol',...
        'BackgroundColor',ui(n).bgc,...
        'Position',pos,...
        'Callback',ui(n).callback);
    else
        uicontrol(hgui,...
        'tag',ui(n).tag,...
        'Style',ui(n).style,...
        'String',ui(n).str,...
        'BackgroundColor',ui(n).bgc,...
        'Position',pos,...
        'Callback',ui(n).callback,...
        'Value',ui(n).value);
    end
    htext_status = findobj('tag','htext_status');
    relax = findobj('tag','relax');
    set(htext_status,'FontSize',18)
    set(relax,'FontSize',16)
    set(hgui,'CloseRequestFcn',@user_closereq)
end

function user_closereq(hgui,~)
% [user_closereq] Brief description:
%   Close request function, runs when GUI is closed.
%   questdlg options: 'Yes', 'Yes and Save Data', 'No/Export'. 
%   All options will send the current data to the base workspace.

data = getappdata(hgui); % gets the current data
% output data structure to workspace
selection = questdlg('QUIT?',...
  'Close Request Function',...
  'Yes','Yes and Save Data','No/Export','Yes'); 
switch selection, 
  case 'Yes', %%%%%%%%%%%%%%%%%%%%%%%% output data structure to workspace
      assignin('base','data',data) % output data
      delete(gcf) % close GUI
  case 'Yes and Save Data' %%%%%%%%%%% save data and output to workspace
      save([data.Name data.Time '.mat'],'data') % save data
      assignin('base','data',data) % output data
      delete(gcf) % close GUI
  case 'No/Export' %%%%%%%%%%%%%%%%%%% cancel but export
      assignin('base','data',data) % output data
  return 
end








