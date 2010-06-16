function [header,data] = parsefile(filename)

% Read in File
id = fopen(filename);
 
% Read in Header Content
readin = fgetl(id);  % get next line
header=regexp(readin, '([^ \t][^\t]*)', 'match');  % split on '\t' into cell array

% Read in Data
data=[];
while ~feof(id) % go until EOF
    i=1;
    data = [data; str2num(fgetl(id))];
    i=i+1;
end
fclose('all');