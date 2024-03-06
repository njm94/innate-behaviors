function show_mov(data,delay,clims,cmap)
if nargin < 2 || isempty(delay), delay = 0.05; end
if nargin < 3 || isempty(clims), clims = [min(data(:)) max(data(:))]; end
if nargin < 4, cmap = 'jet'; end
N = size(data,3);

i = 1;
h = figure('Position',[100 100 1500 800]);
while i < N+1
    if ~ishandle(h)
        break
    end
    imagesc(data(:,:,i)), colormap(cmap); colorbar; caxis(clims)
    title([num2str(i),'/',num2str(N)])
    pause(delay)
    i = i+1;
    if i == N+1
        i = 1;
    end
end