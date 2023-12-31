# ODT-SAS
ODT-SAS is a new semi-automatic method that uses conventional approaches (not based on machine learning) to process optical diffraction tomography (ODT) images of biological cells. It enables 3D segmentation of whole biological cells and 2 types of their organelles: nucleolus and lipid structures.

## INSTALLATION
### REQUIREMENTS
**Base Requirements**
- Matlab R2020b or later
- Graphic card that utilizes CUDA technology, with a Compute Capability >= 2.0
- Installed CUDA driver

**Other Requirements**
- ASTRA Toolbox (Version 1.8 or later: [MATLAB toolbox](https://www.astra-toolbox.com/))
- Spot: A Linear-Operator Toolbox (Version 1.2 or later: [Spot Toolbox](http://www.cs.ubc.ca/labs/scl/spot))

### INSTRUCTION
1. Download the repository
2. Download MATLAB Toolboxes indicated in the requirements and move them to the "ASTRA Tool" folder
3. Run ODT_SAS.m file

_Before running ODT_SAS.m, please read **Code instruction**_

### CITATION
If you find the work useful, please cite the following paper:

M. Mazur and W. Krauze, ["Volumetric segmentation of biological cells and subcellular structures for optical diffraction tomography images"](https://doi.org/10.1364/BOE.498275), Biomed. Opt. Express  14, 5022-5035 (2023).

