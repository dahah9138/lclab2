% File
fpath = "../data/mat/";
f = "Q1_heliknoton_25npp_3dop.mat";

fullpath = strcat(fpath,f);

load(fullpath);

% Convert nn to a list format (x,y,z,nx,ny,nz)
[~,f_strip,e]=fileparts(fullpath);
%fid = fopen(strcat(f_strip,".dat"), 'wt');

[I,J,K,~] = size(nn)

for i=1:I
for j=1:J
for k=1:K
    x = ((i - 1)/(I - 1) - 0.5) * 3;
    y = ((j - 1)/(J - 1) - 0.5) * 3;
    z = ((k - 1)/(K - 1) - 0.5) * 3;
    nx = nn(i,j,k,1);
    ny = nn(i,j,k,2);
    nz = nn(i,j,k,3);
    %fprintf( fid, '%f,%f,%f,%f,%f,%f\n', x, y, z, nx, ny, nz);
end
end
end

%fclose(fid);