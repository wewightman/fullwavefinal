function map = readmap(filename, nx, nz)
%READMAP read the material map
fid = fopen(filename);
vec = fread(fid, inf, 'int8');
fclose(fid);

map = reshape(vec, nz, nx)';
end

