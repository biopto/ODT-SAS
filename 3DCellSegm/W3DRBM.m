function [mask, cellCount] = W3DRBM(REC, n_imm, dx, thresh, option, strel_size, noiseVolume, noiseVolume2, thinFragments, cellThresh) 
% INPUT:
%                  REC: Tomographic reconstruction
%                n_imm: Refractive index of object immersion medium
%                   dx: Object space sample size in XY [um]
%               thresh: Preprocessing threshold
%               option: 1 - remove noise less than noiseVolume [um^3] - default: 500 [um^3]
%                       2 - remove noise less than noiseVolume and greater than noiseVolume2 ([um^3])                              
%           strel_size: Size of structrual element [x y z] for
%                       morphological filtration
%                       ([2 2 2] or [3 3 3] - default: [2 2 2]);
%        thinFragments: Presence of thin fragments (e.g. pseudopods) (0 or 1)
%           cellThresh: Minimum volume of cell [um^3]
%
% OUTPUT:
%                 mask: Mask of cells
%
%% Preprocessing

preREC = preprocessing(REC, n_imm, dx, thresh, option, strel_size, noiseVolume, noiseVolume2);

%% Background and foreground marker
%Background markers
%Binarization of preprocessed tomographic reconstruction
binaryREC = preREC;
binaryREC = imbinarize(binaryREC, 0); %n_imm = 0 after scaling preREC to the interval [0, 1].   

%Morphological opening and dilation
binaryREC = imopen(binaryREC, ones(3,3,3));
binaryREC = imdilate(binaryREC, strel('sphere',5));
binaryREC = imopen(binaryREC, ones(5,5,5));

%Filling holes
for i = 1:size(binaryREC,3)
    binaryREC(:,:,i) = imfill(binaryREC(:,:,i), 'holes');
end

backgroundMarker = ~binaryREC;

%Foreground markers
d = regionprops3(binaryREC,'PrincipalAxisLength');
d = mean(d.PrincipalAxisLength,2);
r = round(max((d/2)));
regionalMaxima = seedPoints(REC, backgroundMarker, r, thinFragments);

%% WATERSHED 3D

%3D Gaussian filtering 
RECSmooth = imcomplement(imgaussfilt3(REC, 4)); 

%Superimpose markers 
imageMin = imimposemin((RECSmooth), backgroundMarker | regionalMaxima, 26);

%Calculate watershed transform 
wsREC = watershed(imageMin);

%% Postprocessing

background = 1; %background marker

%Volume of identified objects
volumeCell = regionprops3(double(wsREC),'Volume'); 
volumeCell.Volume = (dx^3)*volumeCell.Volume; %volume in [um^3]

cellFragm = find(volumeCell.Volume < cellThresh); %fragments of biological cells
A = find(volumeCell.Volume > cellThresh); %biological cells
A(A == background) = []; %remove background marker
fragmA = {};

%Objects identified as a biological cells
if ~isempty(A)
    for n = 1:size(A,1)
        fragm = ismember(wsREC,A(n));
        fragmA{end+1,1} = fragm; 
    end 
end


%Find neighbors of an object
if ~isempty(cellFragm) && ((length(cellFragm) >= 1 && length(A) >= 1) || (length(cellFragm) > 1))
   
    temp = wsREC;
    temp = imdilate(temp,strel('disk',4)); 
    overlap = wsREC & temp;
    pairs = [wsREC(overlap) temp(overlap)]; 
    pairs(any(diff(pairs,[],2) == 0,2),:) = []; %remove duplicate labels
    pairs(pairs(:,1) == background,:) = []; %remove pairs background-object
    pairs = double(pairs);
    pairs = unique(double(pairs), 'rows','stable');
    
    %Sort fragments of biological cell in descending order
   
    for i = 1:size(cellFragm,1)
        cellFragm(i,2) = volumeCell.Volume(double(cellFragm(i,1)));
    end
    cellFragm = sortrows(cellFragm,2,'descend');
    cellFragm(:,2)=[];
    
    fragmPairs = pairs(:,1:2);
    fragmPairs(~ismember(fragmPairs,cellFragm)) = 0;
    fragmPairs = fragmPairs(all(fragmPairs,2),:);
    clear index
    
    for i = 1:size(pairs,1)
        fragmPairs(i,3) = volumeCell.Volume(double(fragmPairs(i,1)));
    end
       
    fragmPairs = sortrows(fragmPairs,3,'descend');
    fragmPairs(:,3) = [];
    
    n = 1;
    n2 = 0;
    while n2 ~= size(fragmPairs,1)
        fragm1 = fragmPairs(n,1);
        n2 = n2 + size(find(fragmPairs(:,1)==fragm1),1); 
        sortPairs = ismember(cellFragm(:,1),fragmPairs(n:n2,2)).*cellFragm(:,1);
        sortPairs(sortPairs==0) = [];
        fragmPairs(n:n2,2) = sortPairs;
        n = n2+1;
    end
    clear fragm1   

    %MERGE FRAGMENTS OF BIOLOGICAL CELL

    while ~isempty(fragmPairs)

        label1 = fragmPairs(1,1); 
        label2 = fragmPairs(1,2); 
        index = find([any(fragmPairs == label2,2)]);

        tempPairs = fragmPairs(index,:);

        column1 = tempPairs(tempPairs ~= label2); 
        column2 = ones(size(tempPairs,1),1).*double(label2); 
        tempPairs = [column1 column2 index];
        overlapFragm = {};
        
        %Determine the common part of adjacent fragments
        for j = 1:size(tempPairs,1)

            fragm1 = ismember(wsREC,tempPairs(j,1)); 
            fragm1 = imdilate(fragm1,strel('disk',4));

            fragm2 = ismember(wsREC,tempPairs(j,2));
            fragm2 = imdilate(fragm2,strel('disk',4));

            overlapFragm{j} = fragm1 & fragm2;
            volumeParts = nnz(overlapFragm{j})*(dx^3); 
            tempPairs(j,4) = volumeParts;

        end
        %find the object that shares the largest volume with another object
        index2 = find([tempPairs(:,4)] == max([tempPairs(:,4)])); 

        if numel(index2) > 1
           index2 = index2(1);
        end

        label1 = tempPairs(index2,1); 
        fragmPairs(ismember(fragmPairs,label2)) = label1;
        pairs(ismember(pairs,tempPairs(index2,2))) = label1;
        fragm = ismember(wsREC,label1) | ismember(wsREC,label2) | overlapFragm{index2};
        volumeFragm = nnz(fragm)*(dx^3);
        wsREC(fragm) = label1;
        
        %Fragments of cell or biological cell
        if volumeFragm < cellThresh
           fragmPairs(tempPairs(index2,3),:) = [];

        else
           fragmPairs(ismember(fragmPairs,label1)) = 0;
           fragmPairs = fragmPairs(all(fragmPairs,2),:);
           A(end+1) = label1;
           fragmA{end+1,1} = fragm;
        end

        fragmPairs(any(diff(fragmPairs,[],2) == 0,2),:) = [];

    end
    
    clear index;
    tempA = num2cell(A);
    label2 = [];

    if ~isempty(A)
        for m = 1:size(A,1)
            index = find([any(ismember(pairs(:,1:2),A(m)),2)]);
            tempA{m,2} = unique([pairs(index,1) pairs(index,2)]);
            tempA{m,2}(ismember(tempA{m,2},A)) = [];
            label2 = [label2 tempA{m,2}']; 
        end

        clear A i
        A = tempA;
        label2 = unique(label2,'stable');
        overlapFragm = {};


        while ~isempty(label2) 
                counter = 0; 
                
                fragm2 = ismember(wsREC,label2(1));
                fragm2 = imdilate(fragm2,strel('disk',4));

                for j = 1:size(A,1)

                    if nnz(ismember([A{j,2}],label2(1))) ~= 0

                        counter = counter+1;
                        label1 = A{j,1}; %label for biological cell
                        fragm1 = fragmA{j}; 
                        fragm1 = imdilate(fragm1,strel('disk',4));
                        overlapFragm{counter,1} = label1;
                        overlapFragm{counter,2} = label2(1); %label for fragment of cell
                        overlapFragm{counter,3} = fragm1 & fragm2; 
                        volumeParts = nnz(overlapFragm{counter,3})*(dx^3);
                        overlapFragm{counter,4} = volumeParts;

                    end

                end
                
                index2 = find([overlapFragm{:,4}] == max([overlapFragm{:,4}]));
                if numel(index2) > 1
                    index2 = index2(1);
                end
                
                label1 = overlapFragm{index2,1}; 
                fragm = ismember(wsREC,label1) | ismember(wsREC,label2(1)) | overlapFragm{index2,3};
                wsREC(fragm) = label1;
                label2(1) = [];
                overlapFragm = {};


        end
    end
end

%% Mask creation
n = 1;
for i = 2:1:max(wsREC(:))
    n = n+1;
    wsREC(ismember(wsREC,i)) = n;  
end

object = 2:1:max(wsREC(:));
mask = ismember(wsREC,object);
mask = imdilate(mask, strel('sphere',4));

for i = 1:size(mask,3)
    mask(:,:,i) = imfill(mask(:,:,i),'holes');
end
mask = imerode(mask, strel('sphere',4));
mask = imerode(mask, strel('disk',1));

cellCount = numel(object);

end




