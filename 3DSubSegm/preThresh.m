function preSegm = preThresh(ROI, n)

    ROI = ROI.^n;
    newROI = rescale(ROI);

    I = newROI; m = size(newROI);
    I = reshape(I, m(1)*m(2)*m(3), 1);

    [counts, binLocations] = imhist(I,1000);
    idx_back = find(counts == max(counts(:))); %Remove background from analysis
    counts(1:idx_back) = [];
    binLocations(1:idx_back) = [];

    best_idx = RosinThreshold(counts);
    thresh = binLocations(best_idx);
    
    preSegm = newROI;
    if thresh > 0.98
        thresh = thresh - 0.1;
    end

    preSegm(preSegm<thresh) = 0;


end