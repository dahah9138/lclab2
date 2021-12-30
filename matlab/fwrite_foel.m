function fwrite_foel(file_name, data, ver)

%% Read header objects again
[Objects, fileIDt, ~] = fread_lmt(file_name);
fclose(fileIDt);
%%

fileID = fopen(file_name, 'wb');
fwrite(fileID, ver, 'int32');
nObj = size(data, 1);
fwrite(fileID, nObj, 'int64');

% Write header objects
for i=1:nObj
    fwrite(fileID, Objects{i, 1}, 'int64');
    fwrite(fileID, convertStringsToChars(Objects{i, 2}));
    fwrite(fileID, Objects{i, 3}, 'int64');
    fwrite(fileID, Objects{i, 4}, 'int64');
end

% Determine type
sz = cell2mat(Objects(1, 3));
if sz == 4
   type = 'float';
else
    type = 'double';
end

% Write data
fwrite(fileID, data{1, 2}, 'int64');
fwrite(fileID, data{2, 2}, 'schar');
fwrite(fileID, data{3, 2}, 'schar');
fwrite(fileID, data{4, 2}, 'int64');
fwrite(fileID, data{5, 2}, 'int');
fwrite(fileID, data{6, 2}, 'schar');
fwrite(fileID, data{7, 2}, type);
fwrite(fileID, data{8, 2}, type);
fwrite(fileID, data{9, 2}, type);
fwrite(fileID, data{10,2}, type);
fwrite(fileID, data{11,2}, type);

fclose(fileID);

end