function [seeds] = seedPoints(REC, initalBG, r, thinFragments)
%INPUT:
%              REC: Tomographic reconstruction
%        initialBG: Initial background of tomographic reconstruction
%                r: Estimated radius value of biological cell
%    thinFragments: Presence of thin fragments (e.g. pseudopods) (0 or 1)
%
%OUTPUT:
%            seeds: Seed points for Watershed
%     
%

%3D Gaussian filtering            
preRECSmooth = imgaussfilt3(REC, 4);

preRECSmooth(initalBG) = 0;

%Opening by reconstruction
if thinFragments == 1
    if r >= 100
        se = strel('disk',round(r/(round(r/100)*15))); 
    else
        se = strel('disk',round(r/15)); 
    end
else
    se = strel('disk',round(r/8)); 
end

marker = imerode(preRECSmooth, se);
morphRecon = imreconstruct(marker,preRECSmooth);

%Regional maxima 
seeds = imregionalmax(morphRecon);

end

