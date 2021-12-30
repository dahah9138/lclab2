function [Data, version] = fread_foel(file_name)


[Objects, fileID, version] = fread_lmt(file_name);

Data = cell(11, 2);

% Read ALL the data

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
tdata = fread(fileID, sz, 'schar');
Data(2,:) = { Objects(2, 2), tdata };

% Relax kind
sz = cell2mat(Objects(3, 3));
tdata = fread(fileID, sz, 'schar');
Data(3,:) = { Objects(3, 2), tdata };

% Iterations
tdata = fread(fileID, 1, 'int64');
Data(4, :) = { Objects(4, 2), tdata };

% voxels
tdata = fread(fileID, 3, 'int');
Data(5, :) = { Objects(5, 2), tdata };
vox = tdata;
N = vox(1) * vox(2) * vox(3);

% Boundaries
tdata = fread(fileID, 3, 'schar');
Data(6,:) = { Objects(6, 2), tdata };


% Cell dims
tdata = fread(fileID, 3, type);
Data(7, :) = { Objects(7, 2), tdata };

% Chirality
tdata = fread(fileID, 1, type);
Data(8, :) = { Objects(7, 2), tdata };

% Relax rate
tdata = fread(fileID, 1, type);
Data(9, :) = { Objects(8, 2), tdata };

% Directors

tdata = fread(fileID, 3 * N, type);
tdata = reshape(tdata, [vox(1) vox(2) vox(3) 3]);
Data(10,:) = { Objects(10, 2), tdata };

% Voltage

tdata = fread(fileID, N, type);
tdata = reshape(tdata, [vox(1) vox(2) vox(3)]);
Data(11,:) = { Objects(11, 2), tdata };

fclose(fileID);

end