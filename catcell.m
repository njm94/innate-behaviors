function newcell = catcell(varargin)
DIM = varargin{1};
assert(isscalar(DIM),'First input must be scalar indicating dimension to concatenate along.')
newcell = [];

for i = 2:nargin
    if iscell(varargin{i})
        for j = 1:length(varargin{i})
            newcell = cat(DIM,newcell,varargin{i}{j});
        end
    else
        newcell = cat(DIM,newcell,varargin{i});
    end
end