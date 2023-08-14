clear
clc
clf
% Open the splay bend file
fname = 'D:/dev/lclab2/data/splay_bend/test.bin';
Ssb_target = 1000;
epsilon = 600;
preimg_alpha = 1.;

fid = fopen(fname, 'r');
vox = fread(fid, [1 3], 'int32');
m = vox(1);
n = vox(2);
p = vox(3);
vol = m*n*p;

Ssb = fread(fid, [1,vol], 'float32');
Ssb = double(reshape(Ssb, [m,n,p]));

diffmag_pos = abs(Ssb - Ssb_target);

[xx,yy,zz] = ndgrid(1:m,1:n,1:p);
fv_pos = isosurface(xx,yy,zz,diffmag_pos,epsilon);

diffmag_neg = abs(Ssb + Ssb_target);

[xx,yy,zz] = ndgrid(1:m,1:n,1:p);
fv_neg = isosurface(xx,yy,zz,diffmag_neg,epsilon);

patch(fv_pos,...
        'FaceColor',[1,1,0],...
        'FaceLighting','gouraud',...
        'FaceAlpha',preimg_alpha,...
        'EdgeColor','none',...
        'SpecularStrength',0.4,...
        'AmbientStrength',0.75)

hold on

patch(fv_neg,...
        'FaceColor',[0,0,1],...
        'FaceLighting','gouraud',...
        'FaceAlpha',preimg_alpha,...
        'EdgeColor','none',...
        'SpecularStrength',0.4,...
        'AmbientStrength',0.75)

view(3)