function [Path, REC, rayXY, n_imm, dx] = loadFiles()

[REC,Path] = uigetfile('*.mat','Select tomographic reconstruction');
fileName=fullfile(Path,REC);
load(fileName);  
clear filename 

if exist('dxo','var') 
    dx = dxo; clear dxo; 
end

if exist('n_immersion','var') 
    n_imm = n_immersion; clear n_immersion; 
end

if exist('RECON','var')
    REC=RECON; clear RECON;
end

if exist('sino_params','var')
    rayXY = zeros(2, size(sino_params(1,:),2));
    rayXY(1,:) = sino_params(3,:)./n_imm .* sino_params(1,:);
    rayXY(2,:) = sino_params(3,:)./n_imm .* sino_params(2,:);
elseif ~exist('rayXY', 'var')   
    rayXY = [];
end

end