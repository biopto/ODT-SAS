function maskSub = segmSubs(REC, maskCell, cellCount, x_tv, n_imm, dx, LS)

%% Seed Points detection
%% AUTOMATIC DETECTION
temp_MC = zeros(size(maskCell{1}));
for i=1:size(mask,1)
    temp_MC = temp_MC + maskCell{i};
end
maskCell = temp_MC;
masked = x_tv;
masked(~maskCell) = 0;

[seedAut] = seedPointSubs(x_tv, masked, cellCount, n_imm, dx);

seedPointAut = regionprops3(seedAut, x_tv, "WeightedCentroid" );
seedPointAut = round(seedPointAut.WeightedCentroid);

%% MANUAL DETECTION

f = manSelection(REC, seedPointAut, [1.33 1.37]); colormap gray;
uiwait(f);

load('seed.mat');

seedPoints = round(seed{1});
dSM = seed{2};                 %maximum nucleolus diameter
nuclRadius = round(dSM/2);     %maximum nucleolus radius

seedPointMan = zeros(size(x_tv));

for i = 1:size(seedPoints, 1)
    seedPointMan(seedPoints(i,2),seedPoints(i,1),seedPoints(i,3)) = 1;
end
seedPointMan = logical(seedPointMan);

%Correct manual seed points selection
props = []; i = 1;
while size(props,1) ~= nnz(seedPointMan)
    seedMan = imdilate(seedPointMan, strel('sphere',nuclRadius-i));       
    props = regionprops3(seedMan, x_tv, "WeightedCentroid" );
    props = round(props.WeightedCentroid);    
    i = i+1;        
end

seedPointMan = zeros(size(x_tv));
for n = 1:size(props,1)
    seedPointMan(props(n,2),props(n,1),props(n,3)) = 1;
end

seedPointMan = logical(seedPointMan);
clear seed seedPoints

%% AUTOMATIC AND MANUAL SEED POINTS

dSA = regionprops3(seedAut,'PrincipalAxisLength');
dSA = mean(dSA.PrincipalAxisLength,2);
rSA = ceil(mean((dSA/2))); %nucleolus radius for automatic detection

param = round(abs(nuclRadius-rSA));
if param == 0
   param = 1;
end
    
%Increase the size of the automatic identified seed points
seedAut_temp = imdilate(seedAut, strel('sphere',param));
labeled = bwlabeln(seedAut_temp, 26);
i = 1;
param2 = 0;
while max(labeled(:)) ~= size(seedPointAut,1)
      param2 = param - i;
      if param2 == 0
         param2 = 1;
         seedAut_temp = imdilate(seedAut, strel('sphere',param2));        
         labeled = bwlabeln(seedAut_temp, 26);
         break;
      end
      seedAut_temp = imdilate(seedAut, strel('sphere',param2));        
      labeled = bwlabeln(seedAut_temp, 26);
      i=i+1;        
end
    
%Select seed points common to automatic and manual detection
seedAM = labeled.*seedPointMan;
idx = seedAM((seedAM(:)>0));
seedROI = ismember(labeled,idx); 
clear labeled

seedPoint = regionprops3(seedROI, x_tv, "WeightedCentroid" );
seedPoint = round(seedPoint.WeightedCentroid); 
if ~isempty(seedPoint)
    seedPoint_temp = zeros(size(x_tv));
    
    for n = 1:size(props,1)
        seedPoint_temp(seedPoint(n,2),seedPoint(n,1),seedPoint(n,3)) = 1;
    end
    
    seedPoint = seedPoint_temp;
else
    seedPoint = seedPointMan;
end

if param2 ~= 0
   seedROI = imdilate(seedROI, strel('sphere',abs(param-param2)));
end
    
seedPointMan(seedROI) = 0;
    
%Add the missing seed points
if nnz(seedPointMan) ~= 0
        
   seedMan = imdilate(seedPointMan, strel('sphere', nuclRadius-1));

   seedMan_temp = imdilate(seedMan, strel('disk',3));
   searchArea = x_tv;
   searchArea(~seedMan_temp) = n_imm;
   seedMan = preThresh(searchArea, 5);

   seedROI = seedROI | seedMan;
           
end
seed = seedPoint | seedPointMan;

if nnz(seedROI.*seed) ~= nnz(seed)
   missed_seed = seed.*(~seedROI);
   missed_seed = logical(missed_seed);
   missed_seed = imdilate(missed_seed, strel('sphere',nuclRadius-1));
   seedROI = seedROI | missed_seed;
end
         

%% ROI of NUCLEOI
%Define the initial ROI of nucleoli
d = regionprops3(seedROI,'PrincipalAxisLength','VoxelList');
d.PrincipalAxisLength = mean(d.PrincipalAxisLength,2);
r = round(d.PrincipalAxisLength/2);
r = min(r(:));
if r < nuclRadius
   param = round(nuclRadius-r) + 5;
else
    param = 5;
end
seedROI = logical(seedROI);
initNucROI = imdilate(seedROI,strel('sphere',param));

temp = x_tv;
temp(~initNucROI) = n_imm;
initNucROI = temp; clear temp;

%Presegmentation - RosinThreshold (histogram)
preNucl = preThresh(initNucROI, 5);

labeled = bwlabeln(preNucl, 26);
tmp = labeled.*seed;
idx = tmp((tmp(:)>0));

preNucl = ismember(labeled,idx); clear labeled idx 
preNucl(preNucl>0)=1;

counter = nnz(tmp);
tmp2 = tmp;
newNucl = zeros(size(x_tv));

%Verify the correctness of segmentation
while counter ~= nnz(seed)
    
    seed_temp = seed;
    seed_temp(logical(tmp)) = 0; 

    labeled = bwlabeln(seedROI,26);
    tmp = labeled.*seed_temp; 
    idx = tmp((tmp(:)>0));
   
    missedNucl = ismember(labeled,idx);
    missedNucl(preNucl) = 0;
    missedNucl = imdilate(missedNucl,strel('sphere',2));

    temp = x_tv;
    temp(~(missedNucl & ~preNucl)) = n_imm;
    missedNucl = temp; clear temp;
    missedNucl = preThresh(missedNucl, 5);
   
    newNucl = newNucl | logical(missedNucl);
    labeled = bwlabeln(newNucl, 26);
    tmp = labeled.*seed_temp; 
    counter = counter + nnz(tmp);
    tmp2 = tmp2 | tmp;

    if counter == 0
        missedNucl = imdilate(logical(missedNucl), strel('sphere',1));
        seedROI(missedNucl) = 0;
        for n = 1:size(seedROI,3)
            seedROI(:,:,n) = imfill(seedROI(:,:,n), 'holes');
        end

    end      
    
end

finalROI = preNucl | (newNucl & ~preNucl);

d = regionprops3(finalROI,'PrincipalAxisLength');
d = mean(d.PrincipalAxisLength,2);
r = (round(d/2));

r = max(r);

if r >= nuclRadius
    se = 2;
else
    se = floor((nuclRadius - r))+1;
end

finalROI = imdilate(finalROI, strel('disk',se));
finalROI=~finalROI & imdilate(finalROI, strel('sphere',1));
x_tv = rescale(x_tv);

%3D watershed - final segmentation of nucleoli
imageMin = imimposemin(imcomplement(x_tv), finalROI | seed |  ~masked, 26);
wsREC = watershed(imageMin);

tmp = double(wsREC).*seed; 
idx = tmp((tmp(:)>0));
maskNucl = ismember(wsREC, idx);

%% Lipid structures segmentation 
if LS == 1
    %Remove background and nucleoli from analysis
    maskedREC = REC;
    maskedREC(~maskCell) = 0;
    limitedArea = maskedREC;
    limitedArea(limitedArea == 0) = n_imm;
    limitedArea(maskNucl) = n_imm;

    limitedArea = rescale(limitedArea); %interval [0 1]

    %Enhance contrast
    se = strel('sphere',3);
    J = imsubtract(imadd(limitedArea,imtophat(limitedArea,se)),imbothat(limitedArea,se));

    %Initial ROI of LS
    preLipid = logical(preThresh(J, 5));
    preLipid = imdilate(preLipid, strel('sphere',3));
    J(~preLipid) = 0;

    %Final segmentation
    maskLipid = logical(preThresh(J, 2));
    maskLipid = imdilate(maskLipid,strel('disk',1));
else
    maskLipid = [];
end
    

%% Nucleoli and LS masks
varNames = {'Nucleoli','LS'};
varTypes = {'logical','logical'};

maskSub = table('Size',[1 2],'VariableTypes',varTypes,'VariableNames',varNames);
maskSub.Nucleoli = {maskNucl};
maskSub.LS = {maskLipid};

end