function [x,y] = get_neuron_location_in_dorsal_space(x0,y0,roi_tform, linear_tform, wfi_tform)   
    
    [x1, y1] = transformPointsForward(roi_tform.tform, x0, y0);        
    [x2, y2] = transformPointsForward(linear_tform.tform, x1, y1);        
    [x, y] = transformPointsForward(wfi_tform.tform, x2, y2);


end