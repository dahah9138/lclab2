function [data, num_points, interaction_type] = load_interaction_data2(file)

% don't touch
s = dir(file);
double_sz = 8;
datastruct_size = 10 * double_sz;
num_points = (s.bytes - 1) / datastruct_size;

fID = fopen(file, 'r');
interaction_type = fread(fID, 1, 'char*1');
data = fread(fID, [10 num_points], 'double');
fclose(fID);

end