function [seed] = seedPointSubsv5(x_tv, masked, cellCount, n_imm, dx)

%Thresholding
RI_2D = zeros(size(masked,3),2);
 
for i = 1:size(masked,3)
    masked_2D = masked(:,:,i);
    RI_2D(i,1) = i;
    RI_2D(i,2)  =sum(masked_2D(:))/nnz(masked_2D(:));
end
x_tv(~masked) = n_imm;
x_tv(x_tv <= max(RI_2D(:,2))) = n_imm;

%Estimation of the nucleolus radius 
cytoplasm_volume = nnz(masked)*(dx^3); %[um^3]
cytoplasm_volume = cytoplasm_volume/cellCount;
percent = (cytoplasm_volume/((numel(x_tv))*(dx^3)))*100; 
cytoplasm_radius = ((cytoplasm_volume/pi)*(3/4))^(1/3); 
nucleolus_radius = (0.9/5.8)*cytoplasm_radius; %[um]
nucleolus_radius = round((nucleolus_radius/dx)); 

%Morphological reconstruction
if percent < 1
   m=3;
else
   m=2;
end
    
if round(nucleolus_radius/m) <= 1
   se = strel('disk',2); 
else
   se = strel('disk',round(nucleolus_radius/m)); 
end

x_tv = rescale(x_tv);
marker = imerode(x_tv,se);
morphx_tv = imreconstruct(marker,x_tv);


seed = imregionalmax((morphx_tv));

for i = 1:size(seed,3)
    seed(:,:,i) = imfill(seed(:,:,i),'holes');
end   

end