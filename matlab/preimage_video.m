% Plotter for preimage, xsections
clear 
tic
% figure
clf
hold all
delete(findall(gcf,'Type','light'))

file = '../data/mat/loop_torus_cropped';

% Cone error
eps = 0.15;
fps = 20;
thetaframes = 20;
phiframes = 50;
numframes = thetaframes * phiframes;

% Sampling rate for isosurfaces of preimages
sample = 1; 
preimg_alpha = 1;

video = VideoWriter('../data/videos/test');
video.FrameRate = fps;
open(video);


file = strcat(file, '.mat');

try
    load(file,'nn');
    disp('loaded nn')
catch    
    disp('failed to load nn')
end

% Cut off some of nn

[n0,m0,p0,~] = size(nn);

p01 = int32(p0 *3/9);
p02 = int32(p0 *6/9);

n01 = int32(n0 *2/9);
n02 = int32(n0 *7/9);

m01 = int32(m0 * 2/9);
m02 = int32(m0 * 7/9);

cutnn = nn(n01:n02,m01:m02,1:p02,:);
%nn = cutnn;

viewangle = [0,0];
xsection_alpha = 0.4;

% Sampling rate for the xz-xsection
ratexs = 4; 

% parametrized reference vector
numpreimages = 1;
%thetad = repmat([60, 60],numframes,1);
thetad = zeros(numframes, numpreimages);
%thetad(:,1) = 0;
%thetad(:,2) = linspace(0, 180, numframes);

phid = zeros(numframes, numpreimages);
%phid(:,1) = linspace(0, 360, numframes);
%phid(:,2) = linspace(180, 540, numframes);

% Map id = phi * THETA + theta
for p = 1:phiframes
    ph0 = (p-1) * 2*pi / (phiframes - 1) * 180/pi;
   for t = 1:thetaframes
        id = (t-1) * phiframes + p;
        thetad(id,1) = (t-1) * pi / (thetaframes - 1) * 180/pi;
        phid(id,1) = ph0;
   end
end


epsilon = eps*ones(1,numpreimages);

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

for k=1:numframes
    clf('reset')
    


    phi_cur = phid(k, 1);
    theta_cur = thetad(k, 1);
    
    
    
    % Draw the sphere
    kk = 10;
    map = hsv(2^kk);
    radius = sqrt(ma + na + pa);
    
    
    for cone=1:numpreimages

        % Draw the cone on the sphere
        [x1,y1,z1] = cylinder([1 0],128);

        xc = x1;
        yc = y1;

        zc = radius*(z1+.8);

        % set up rotation matrix1:
        rotationMatrix1 = [cos(theta(k,cone))  0  sin(theta(k,cone));...
                           0                  1  0;...
                           -sin(theta(k,cone)) 0  cos(theta(k,cone))];
        % set up rotation matrix2:
        rotationMatrix2 = [cos(phi(k,cone))  -sin(phi(k,cone))  0;...
                           sin(phi(k,cone))  cos(phi(k,cone))   0;...
                           0                0                 1];
        % get points at the two rings and rotate them separately about x:

        positionNew1 = rotationMatrix2*rotationMatrix1*[xc(1,:); yc(1,:); zc(1,:)];
        positionNew2 = rotationMatrix2*rotationMatrix1*[xc(2,:); yc(2,:); zc(2,:)];

        xc = [positionNew1(1,:); positionNew2(1,:)]+ 1.25*m;%-2*radius;
        yc = [positionNew1(2,:); positionNew2(2,:)]+ 1.25*n;%-2*radius;
        zc = [positionNew1(3,:); positionNew2(3,:)]+ 1.25*p;%-2*radius;

        arrow_color = zeros(1,3);
        arrow_color(1) = (mod(phi(k,cone), 2*pi))/(2*pi); 
        if arrow_color(1)<0
            arrow_color(1) = (mod(phi(k,cone), 2*pi)+pi)/(2*pi); 
        end
        arrow_color(2) = theta(k,cone)/(pi/2);
        arrow_color(3) = 2 - theta(k,cone)/(pi/2);

        arrow_color(arrow_color>1) = 1;
        arrow_color = hsv2rgb(arrow_color);

        
        s = surf(xc,yc,zc,'EdgeColor','none','FaceColor',arrow_color);
        fvc = surf2patch(s);
        delete(s)
        hold on
        patch(fvc,...
            'FaceAlpha',1,...
            'FaceColor',arrow_color,...
            'FaceLighting','flat',...
            'EdgeColor','none',...
            'SpecularStrength',0.8,...
            'AmbientStrength',0.55);
        
        hold on
        
        phixs = sin(theta(k,cone)).*cos(phi(k,cone));
        py = sin(theta(k,cone)).*sin(phi(k,cone));
        pz = cos(theta(k,cone));

        diffmag = sqrt((nn2(:,:,:,1)-phixs).^2+...
            (nn2(:,:,:,2)-py).^2+(nn2(:,:,:,3)-pz).^2);

        disp('Finished diffmag')

        fv = isosurface(xx,yy,zz,diffmag,epsilon(cone));

    disp('Created a isosurface')
        
        if numel(fv.vertices)>3

            color = zeros(1,3);
            color(1) = mod(phi(k,cone), 2*pi)/2/pi;
            color(2) = theta(k,cone)/(pi/2);
            color(3) = 2 - theta(k,cone)/(pi/2);
            color(color > 1) = 1;
            color = hsv2rgb(color);

            patch(fv,...
                'FaceColor',color,...
                'FaceLighting','gouraud',...
                'FaceAlpha',preimg_alpha,...
                'EdgeColor','none',...
                'SpecularStrength',0.4,...
                'AmbientStrength',0.75);

            disp('finished patching')
        %     hold off


        end

        disp('Made a preimage')

    end
    
    %%%%%%

    hold on


    [X,Y,Z]= sphere(500);

    [ms,ns] = size(X);
    [azi,ele,~] = cart2sph(X,Y,Z);
    % Z(Z >= cos(bound_d*pi/180)) = nan; % remove the region with theta<bound

    azi_rs = reshape(azi,[ms*ns,1]);
    azi_rs(azi_rs<0) = azi_rs(azi_rs<0)+2*pi;
    the_rs = pi-(reshape(ele,[ms*ns,1])+pi/2);
    col = zeros(ms*ns,3);

    % hsv: s = 0 at the north pole; v = 0 at the south pole;
    for i = 1:ms*ns
        col(i,1) = (azi_rs(i))/(2*pi); 
        col(i,2) = the_rs(i)/(pi/2);
        col(i,3) = 2 - the_rs(i)/(pi/2);

    end
     % Translate
    X = radius * X + 1.25*m;%-2*radius;
    Y = radius * Y + 1.25*n;%-2*radius;
    Z = radius * Z + 1.25*p;%-2*radius;

    col(col > 1) = 1;
    col_rs = reshape(col,[ms,ns,3]);
    % figure
    s1=surf(X,Y,Z,hsv2rgb(col_rs)); % color-code according to theta
    s1.EdgeColor = 'none';
    s1.SpecularStrength = 0.01;
    s1.AmbientStrength= 0.7;
    s1.FaceAlpha = 1;

    colormap(map)
    %%%%%%

    box
    ax = gca;
    ax.BoxStyle = 'full';
    axis off

    set(gcf,'Color','w')
    text(ma/2,-na/2,sprintf("(theta, phi) = (%.1f, %.1f)", theta_cur, phi_cur));
    ax.XTick = [];
    ax.YTick = [];
    ax.ZTick = [];
    axis equal
    xlim([1 1.5*m])
    ylim([1 1.5*n])
    zlim([1 1.5*p])
    view(viewangle)
    camlight
    camproj('p')
    % %%
    % %% axis arrows
    %hold all
    axis off
    box on
    xlim([-2 1.5*m+1])
    ylim([-2 1.5*n+1])
    zlim([-2 1.5*p+1])


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

    scale = 1.0;
    small_bd = (m+1)/2-(m+1)/2*scale;
    large_bd = (m+1)/2+(m+1)/2*scale;
    xlim([small_bd-2 2*large_bd+1])
    ylim([small_bd-2 2*large_bd+1])
    zlim([-2 2*p+1])
    linewidth = 1;
    
    frame = getframe(gcf);
    writeVideo(video, frame);
end

close(video);
