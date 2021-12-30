function Data = fread_fd(file_name)


[Objects, fileID, ~] = fread_lmt(file_name);

% Now read and format the data
% In this case, extract voxels, cell_dims, iterations, and directors
% To get directors from Data, use cell2mat(Data(4,2))

Data = cell(4, 2);

% Scalar size
sz = cell2mat(Objects(1, 3));
tdata = fread(fileID, 1, 'int64');
Data(1, :) = { Objects(1, 2), tdata };

if sz == 4
   type = 'float';
else
    type = 'double';
end

% LC type
sz = cell2mat(Objects(2, 3));
fread(fileID, sz, 'schar');

% Relax kind
sz = cell2mat(Objects(3, 3));
fread(fileID, sz, 'schar');

% Iterations
tdata = fread(fileID, 1, 'int64');
Data(1, :) = { Objects(4, 2), tdata };

% voxels
tdata = fread(fileID, 3, 'int');
Data(2, :) = { Objects(5, 2), tdata };
vox = tdata;
N = vox(1) * vox(2) * vox(3);

% Boundaries
fread(fileID, 3, 'schar');

% Cell dims
tdata = fread(fileID, 3, type);
Data(3, :) = { Objects(7, 2), tdata };

% Chirality
fread(fileID, 1, type);

% Relax rate
fread(fileID, 1, type);

% Directors

tdata = fread(fileID, 3 * N, type);
tdata = reshape(tdata, [vox(1) vox(2) vox(3) 3]);
Data(4,:) = { Objects(10, 2), tdata };

end