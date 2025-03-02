function [data, num_points, interaction_type] = load_interaction_data2(file)

% don't touch
df = dir(file);
double_sz = 8;
nFields = 4;
datastruct_size = nFields * double_sz;
fbytes = df.bytes;
num_points = (fbytes - 1) / datastruct_size;

fID = fopen(file, 'r');
interaction_type = fread(fID, 1, 'char*1');
data = fread(fID, [nFields num_points], 'double');
fclose(fID);

end