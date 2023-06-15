% Authors: Wojciech Krauze
% Contact: wojciech.krauze@pw.edu.pl
% Affiliation: Warsaw University of Technology, Institute of Micromechanics and Photonics

function x_tv = generate_xtv(REC, n_imm, rayXY, lambda)
% INPUT:
%                  REC: masked tomographic reconstruction
%           rayXY(1,:): x-coordinates of illumination vectors
%           rayXY(2,:): y-coordinates of illumination vectors
%               lambda: weighting factor which controls the "intensity" of
%                       the TV regularization
%
% OUTPUT:
%                 x_tv: regularized tomographic reconstruction (TV
%                       regularization)
%

if ~isempty(rayXY)
    rx = rayXY(1,:);
    ry = rayXY(2,:);
else
    % syntetic generation of projection directions
    xa = linspace(0,2*pi-2*pi/90,90);
    A = 0.7;
    rx = A*sin(xa);
    ry = A*cos(xa);
end

% Generate illumination versors from illumination directions
vectors = rays_to_vectors_new(rx,ry, 'fixed');

% ASTRA parameters
N_CP = 250;
vol_geom = astra_create_vol_geom(N_CP, N_CP, N_CP);
proj_geom = astra_create_proj_geom('parallel3d_vec', N_CP, N_CP, vectors);

data = ndresize(REC-n_imm,[N_CP, N_CP, N_CP],'cubic');
A = opTomo('cuda', proj_geom, vol_geom);

projdata = A*data(:);

TV3D = astra.tv.opTV3D(N_CP);

starting_point = zeros(N_CP,N_CP,N_CP);
scale = 0.05*4;

x_tv = astra.tv.chambolle_pock(...
    scale*A, TV3D, scale*projdata, 0, lambda, true, starting_point(:));
x_tv = single(reshape(x_tv,[N_CP N_CP N_CP]));
x_tv = x_tv+n_imm;

end