function [data, num_points, interaction_type] = load_interaction_data(file)

% don't touch
s = dir(file);
double_sz = 8;
datastruct_size = 5 * double_sz;
num_points = (s.bytes - 1) / datastruct_size;

fID = fopen(file, 'r');
interaction_type = fread(fID, 1, 'char*1');
data = fread(fID, [5 num_points], 'double');
fclose(fID);

end