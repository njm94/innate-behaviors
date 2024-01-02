function visualize_registration(varargin)
% input 1 wfi frame
% input 2 save gif boolean
    fig=figure;
    imagesc(varargin{1}); colormap gray, axis off
    hold on
    if ~isempty(varargin{2})
        drawnow
        frame = getframe(fig);
        im{1} = frame2im(frame);
    end

        for i = 3:length(varargin)
            imAlpha=ones(size(varargin{i}));
            imAlpha(varargin{i}==0)=0;
            imagesc(varargin{i}, 'AlphaData', imAlpha)
        end

    if ~isempty(varargin{2})
        drawnow
        frame = getframe(fig);
        im{2} = frame2im(frame);
    

        filename = varargin{2}; % Specify the output file name
        for idx = 1:2
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.5);
            else
                imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.5);
            end
        end
    end

end