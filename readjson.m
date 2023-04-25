function jsonstruct = readjson(filename)
%READJSON return the contents of a json file as a structure
fid = fopen(filename);
raw = fread(fid, inf)';
str = char(raw);
fclose(fid);

jsonstruct = jsondecode(str);
end

