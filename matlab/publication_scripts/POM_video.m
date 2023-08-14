clc
clf
clear

filepath = 'D:/dev/lclab2/data/knot/test';
fileprefix = '';
fps = 25;
upsample = 10;

a=dir([filepath '/POM*.bmp']);
nFrames = size(a,1);
videoname = 'POM_video';
vid = VideoWriter([filepath '/' videoname]);
vid.FrameRate = fps;
vid.Quality = 100;
open(vid)

for k=1:nFrames
    clf('reset')
    f = gcf;
    f.Position = [100 100 840 800];

    file = [filepath '/' fileprefix 'POM_' num2str(k) '.bmp'];
    A = imread(file);

    % Interpolate image data
    sz = size(A);
    xg = 1:sz(1);
    yg = 1:sz(2);
    F = griddedInterpolant({xg,yg},double(A));
    xq = (0:1/upsample:sz(1))';
    yq = (0:1/upsample:sz(2))';
    vq = uint8(F({xq,yq}));

    vq = imgaussfilt(vq, 1);
    image(vq)

    frame = getframe(f);
    writeVideo(vid, frame);

end

close(vid)