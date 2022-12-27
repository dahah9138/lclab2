function f_preimage(nn, sample, thetad, phid)

% Cone error
eps = 0.2;

preimg_alpha = 1;

epsilon = eps*ones(size(thetad));

theta = thetad/180*pi;
phi=phid/180*pi;

[m,n,p,~] = size(nn);

xx = linspace(1,m,m*sample);
yy = linspace(1,n,n*sample);
zz = linspace(1,p,p*sample);

[xx,yy,zz] = ndgrid(xx,yy,zz);

nn_int(:,:,:,1) = interp3(nn(:,:,:,1),yy,xx,zz,'cubic');
nn_int(:,:,:,2) = interp3(nn(:,:,:,2),yy,xx,zz,'cubic');
nn_int(:,:,:,3) = interp3(nn(:,:,:,3),yy,xx,zz,'cubic');

% % normalize director field after interpolation
modnn = sqrt(nn_int(:,:,:,1).^2+nn_int(:,:,:,2).^2+nn_int(:,:,:,3).^2);
nn_int(:,:,:,1) = nn_int(:,:,:,1)./modnn;
nn_int(:,:,:,2) = nn_int(:,:,:,2)./modnn;
nn_int(:,:,:,3) = nn_int(:,:,:,3)./modnn;

for cone=1:length(theta)
    
    phixs = sin(theta(cone)).*cos(phi(cone));
    py = sin(theta(cone)).*sin(phi(cone));
    pz = cos(theta(cone));
    
    diffmag = sqrt((nn_int(:,:,:,1)-phixs).^2+...
        (nn_int(:,:,:,2)-py).^2+(nn_int(:,:,:,3)-pz).^2);
    
    fv = isosurface(xx,yy,zz,diffmag,epsilon(cone));
    
    fprintf('Isosurface (%d/%d)\n',cone,length(theta));
    
    
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
    

end

end