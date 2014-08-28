function vol_out = range(vol, min_val, max_val)

    
vol_out = zeros(size(vol));

max_vols = max(vol(:));
min_vols = min(vol(:));
    

for i = 1:size(vol,3)
    
    vol_out(:,:,i) = min_val + (vol(:,:,i) - min_vols).*((max_val - min_val)/(max_vols - min_vols));
    
end

% vol_out = vol./max_val;