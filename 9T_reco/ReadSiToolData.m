% Rolf Pohmann %
% Max Planck Institute for Biological Cybernetics %

function data = ReadSiToolData( path )
%Reads images stored by SiTools as .rms files %Rolf Pohmann%
%   

[pathname, filename, ext] = fileparts(path);
[tok,remain] = strtok(filename,'_');
compl = false;
ind = 0;
type = 'none';

while (numel(tok)>0)
    t = str2num(tok);
    if (numel(t)==0)
        if strcmp(tok, 'flreal')
            type = 'float32';
        elseif strcmp(tok, 'flcompl')
            type = 'float32';  
            compl = true;
        end
    else
        ind = ind+1;
        dims(ind) = t;
    end
        
    [tok,remain] = strtok(remain,'_');
end

f = fopen(path);
nel = 1;
for cnt = 1:numel(dims)
    nel = nel*dims(cnt);
end
if compl
    data = fread(f,[2,nel],type);
    %nel = nel*2;
else
    data = fread(f,nel,type);
end
%data = fread(f,nel,type);
fclose(f);
if compl
    %s = numel(data);
    %data = reshape(data,[2,nel/2]);
    data = complex(data(1,:),data(2,:));
end
data = reshape(data,dims);
u = zeros(numel(dims),1);
u(1) = round(dims(1)/2);
u(2) = round(dims(2)/2);
u(3) = round(dims(3)/2);
data = circshift(data,u);
end

