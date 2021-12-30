function Data = fread_RBF(file_name)


[Objects, fileID, ~] = fread_lmt(file_name);

% Now read and format the data
% In this case, extract
% - iterations
% - cell_dims
% - directors
% - positions
% - voltage
% - active node indices

% To get directors from Data, use cell2mat(Data(3,2))
% To get positions from Data, use cell2mat(Data(4,2))
% To get voltage from Data, use cell2mat(Data(5,2))
% To get active_nodes from Data, use cell2mat(Data(6,2))

Data = cell(6, 2);

sz = fread(fileID, 1, 'int64');

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

% Boundaries
fread(fileID, 3, 'schar');

% Cell dims
tdata = fread(fileID, 3, type);
Data(2, :) = { Objects(6, 2), tdata };

% Chirality
fread(fileID, 1, type);

% Relax rate
fread(fileID, 1, type);

% Num nodes
N = fread(fileID, 1, 'int64');

% Num active nodes
Nact = fread(fileID, 1, 'int64');

% Num nearest neighbors (per node)
fread(fileID, 1, 'int64');

% Directors
tdata = fread(fileID, 3 * N, type);
tdata = reshape(tdata, [N 3]);
Data(3,:) = { Objects(12, 2), tdata };

% Positions
tdata = fread(fileID, 3 * N, type);
tdata = reshape(tdata, [N 3]);
Data(4,:) = { Objects(13, 2), tdata };

% Voltage
tdata = fread(fileID, N, type);
Data(5,:) = { Objects(14, 2), tdata };

% Active nodes
tdata = fread(fileID, Nact, 'int64');
Data(6,:) = { Objects(15, 2), tdata };


end