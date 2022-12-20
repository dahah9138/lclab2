function fwrite_mat_to_foel(file_in, file_out)

% Get nn from input file
load(file_in, 'nn');

version = 1;

sz_int = 4; % bytes
sz_int64 = 8; % bytes
sz_char = 1; % bytes
sz_double = 8; % bytes
type = 'double'; % scalar type

if type == 'double'
    sz_scalar = 8;
else
    sz_scalar = 4;
end

voxels = size(nn(:,:,:,1));
nn_vol = 3*voxels(1) * voxels(2) * voxels(3);
nObj = 3;
Objects = cell(nObj,4);

Objects(1,:) = {6*sz_char,"Voxels",3*sz_int,0};
Objects(2,:) = {9*sz_char,"Directors",nn_vol*sz_double,1};
Objects(3,:) = {11*sz_char,"Scalar size",sz_int64,2};

fileID = fopen(file_out, 'wb');

% Beginning file information
% Header version (Only one version V = 1)
fwrite(fileID, version, 'int32');
% Number of objects to read from file
fwrite(fileID, nObj, 'int64');

% Write header objects
for i=1:nObj
    fwrite(fileID, Objects{i, 1}, 'int64');
    fwrite(fileID, convertStringsToChars(Objects{i, 2}));
    fwrite(fileID, Objects{i, 3}, 'int64');
    fwrite(fileID, Objects{i, 4}, 'int64');
end

% Write data
fwrite(fileID, voxels, 'int32');
fwrite(fileID, nn, type);
fwrite(fileID, sz_scalar, 'int64');

fclose(fileID);

end