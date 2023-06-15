close all
% clear all
tic

[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(path));

%% Select type of segmentation
% 'Cell';                 %3D segmentation of whole cell
% 'Organelle';            %3D segmentation of nucleolus and lipid droplets
% 'Cell_Organelle';       %3D segmentation - cell + organelle


type = 'Organelle';


%% Load files
[Path, REC, rayXY, n_imm, dx] = loadFiles();

%% Parameters for 3D cell segmentation
%1) Preprocessing:
thresh = 1.342; 
strel_size = [];            %[2 2 2] or [3 3 3]
option = 2;                 %option = 1 or option = 2; 
noiseVolume = 2000;         %[um^3]; 
noiseVolume2 = 10000;       %[um^3] (if option == 2)

%2)Seed points
thinFragments = 1;          %set thinFragments to 1 if the analyzed biological cells 
                            %have complex shape containing thin fragments
                            %such as pseudopodia

%3) Postprocessing:
cellThresh=15000;           %[um^3];

%% Parameters for 3D organelle segmentation
LS = 1;                     %1 - lipid structure segmentation
                            %0 - only nucleolus segmentation
cellCount = 1;              %number of cells analyzed (select only in 'Organelle' type of segmentation)

%% 3D Segmentation

switch type
    case 'Cell' 
       [maskCell, cellCount] = W3DRBM(REC, n_imm, dx, thresh, option, strel_size, noiseVolume, noiseVolume2, thinFragments, cellThresh);
        masked = REC.*maskCell;
        File = strcat(Path,'maskCell.mat');
        save(File,'maskCell', '-v7.3')

    case 'Organelle'
        [mask,Path] = uigetfile('*.mat','Select cell mask');
        fileName = fullfile(Path,mask);
        maskCell = load(fileName);
        
        x_tv = generate_xtv(REC, n_imm, rayXY, 0.004);
        x_tv = ndresize(x_tv, size(REC), 'linear');

        maskSub = segmSubs(REC, maskCell, cellCount, x_tv, n_imm, dx, LS);
        File = strcat(Path,'maskSub.mat');
        save(File,'maskSub', '-v7.3')

    case {'Cell_Organelle'}
       [maskCell, cellCount] = W3DRBM(REC, n_imm, dx, thresh, option, strel_size, noiseVolume, noiseVolume2, thinFragments, cellThresh);
        
        x_tv = generate_xtv(REC, n_imm, rayXY, 0.004);
        x_tv = ndresize(x_tv, size(REC), 'linear');

        maskSub = segmSubs(REC, maskCell, cellCount, x_tv, n_imm, dx, LS);
        File = strcat(Path,'maskCellSub.mat');
        save(File,'maskCell','maskSub', '-v7.3')

    otherwise
        warning('Unexpected type of segmentation.')
end

