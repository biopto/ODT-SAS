function preREC = preprocessing(REC, n_imm, dx, thresh, option, strel_size, noiseVolume, noiseVolume2)
 
if ~exist('thresh','var') || isempty(thresh)
    error('Threshold for preprocessing is required!');
end

%Thresholding:
REC(REC < thresh) = n_imm; 

newREC = rescale(REC); % interval [0,1]
        
%Morphological filtration:   
if ~exist('strel_size','var') || isempty(strel_size) || ~isequal([2 2 2], strel_size) || ~isequal([3 3 3], strel_size)
    strel_size = [2 2 2];
end
    
se = ones(strel_size);
newREC = (imclose(imopen(newREC,se), se));
    
%Labelling:
binaryREC = newREC;
binaryREC(binaryREC > 0) = 1;

for i = 1:size(binaryREC, 3)
    binaryREC(:,:,i) = imfill(binaryREC(:,:,i),'holes');
end

labeledREC = bwlabeln(binaryREC, 26);

volumeRegion = regionprops3(labeledREC, 'Volume'); %[voxels]
voxelCount = volumeRegion.Volume;
volumeRegion.Volume = (dx^3) * voxelCount; %[um^3]
    
if ~exist('noiseVolume','var') || isempty(noiseVolume)
   noiseVolume = 500;
end
    
%Remove artifacts:
if option == 1
   noise = find(volumeRegion.Volume < noiseVolume);
        
elseif option == 2
    noise = find(volumeRegion.Volume < noiseVolume | volumeRegion.Volume > noiseVolume2);
end
    
noisyReg = ismember(labeledREC, noise);
    
preREC = newREC;
preREC(noisyReg) = 0;

end