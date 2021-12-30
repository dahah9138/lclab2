% Author: Darian Hall
% Date: July 23, 2021
% Reads data from .lmt file written by lclab2

function [Objects, fileID, Version] = fread_lmt(file_name)

fileID = fopen(file_name, 'r');

% First read the header
% ---
% Version
Version = fread(fileID, 1, 'int32');

% ---

% ---
% Number of Objects
nObj = fread(fileID, 1, 'int64');
% ---

% ---
% Get objects
Objects = cell(nObj, 4);
for i=1:nObj
    vsize = fread(fileID, 1, 'int64');
    vname = convertCharsToStrings(cast(fread(fileID, vsize, 'schar'), 'char'));
    bytesize = fread(fileID, 1, 'int64');
    loc = fread(fileID, 1, 'int64');
    Objects(i, :) = { vsize, vname, bytesize, loc };
end

end