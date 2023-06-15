%
% best_idx = RosinThreshold(hist_img)
%
% Implementation of Rosin Thresholding. 
% Compute the Rosin threshold for an image
% Takes histogram of an image filtered by any edge detector as as input
% and return the index which corresponds to the threshold in histogram
% 
% REF: "Unimodal thresholding" by Paul L. Rosin (2001)
%
function best_idx = RosinThreshold(hist_img)

% find best threshold

[peak_max, pos_peak] = max(hist_img);
p1 = [pos_peak, peak_max];

% find last non-empty bin
ind_nonZero = find(hist_img>0);
last_zeroBin = ind_nonZero(end);
p2 = [last_zeroBin, hist_img(last_zeroBin)];
best_idx = -1;
max_dist = -1;
for x0 = pos_peak:last_zeroBin
    y0 = hist_img(x0);
    a = p1 - p2;
    b = [x0,y0] - p2;
    cross_ab = a(1)*b(2)-b(1)*a(2);
    d = norm(cross_ab)/norm(a);
    if(d>max_dist)
        best_idx = x0;
        max_dist = d;
    end
end