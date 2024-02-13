function movingReg = evaluate_tform(expt_frame, tform, template_frame)
if nargin<3 || isempty(template_frame), template_frame=[]; end

    Rfixed = imref2d(size(expt_frame(:,:,1)));
    movingReg = imwarp(expt_frame, tform, "OutputView", Rfixed);
    
    if nargout < 1 || ~isempty(template_frame)
        figure,
        subplot(1,2,1)
        imshowpair(template_frame, expt_frame), title('Before registration')
        subplot(1,2,2)
        imshowpair(mouse_ref_image, movingReg), title('After registration')
    end
end