% algorithms only for parallel geometry, and only for 2D CUDA
classdef opRecon < opSpot
    
   properties (Access = private)
       funHandle
       proj_geo
       vol_geo
       proj_handle
       model_type
       sino_id   % to indicate where the sinogram is stored in memory
       rec_id
       initial_vol_id
       alg_id
       % handle for freeing memory
       sino_handle
       rec_handle
%        initial_vol_id_handle
       num_iter
   end
   
   properties (SetAccess=private, GetAccess=public)
       % publice read-only size of Recon operator is necessary for a user to appropriatly reshape the 'x' 
       % matrix into vector
       proj_size
       rec_size
   end
   
   properties (SetAccess=public, GetAccess=public)
       x0
   end
    
   methods (Access = public)
       % CONSTRUCTOR, rec_size is a scalar, nangles is a scalar with number
       % of angles, angle_range is a scalar with angle range in radians
       function op = opRecon(type,rec_size,angle_range_beg,angle_range_end,nangles,numiter)
          if nargin<4
             error('Not enough input parameters: type,nsize,msize,nangles,numiter') 
          end
%           if (nargin<5 && isempty(type))
%              type='CUDA'; 
%           else
%               if (nargin<5) error('Initial parameters are missing')
%               end
%           end
%           if ~strcmp(type,'CUDA')
%              error('For now the only available type is CUDA') 
%           end
          
          proj_geo = astra_create_proj_geom('parallel',1.0,rec_size,linspace2(angle_range_beg,angle_range_end,nangles)) ;
          vol_geo = astra_create_vol_geom(rec_size,rec_size) ;
          
          % image data
          proj_size = astra_geom_size(proj_geo);
          vol_size  = astra_geom_size(vol_geo);
          
          % handle, for freeing data if object is deleted
%           rec_handle = astra.data2d_handle(rec_id2);
          % give me op object created within opSpot superclass with given arguments (that are specified by 
          % opSpot superclass: type,m,n where type is the name of a new Spot operator and m,n is the size of 
          % that operator)
          op = op@opSpot('opRecon',prod(vol_size),prod(proj_size)) ;
          
          sino_id = astra_mex_data2d('create', '-sino', proj_geo, 0);
          rec_id  = astra_mex_data2d('create', '-vol', vol_geo, 0);
%           initial_vol_id  = astra_mex_data2d('create', '-vol', vol_geo, 0);
          
          op.sino_handle = astra.data2d_handle(sino_id);
          op.rec_handle  = astra.data2d_handle(rec_id);
%           op.initial_vol_id_handle = astra.data2d_handle(initial_vol_id);
          
          % handle, for freeing data if object is deleted
          %                 sino_handle = astra.data2d_handle(sino_id2);
          
          % Configuration for astra fp algorithm
          cfg_Recon = astra_struct(type);
          cfg_Recon.ProjectionDataId = sino_id;
%           cfg_Recon.VolumeDataId     = initial_vol_id;
          cfg_Recon.ReconstructionDataId = rec_id;
          cfg_Recon.option.GPUindex  = 0;
          op.alg_id = astra_mex_algorithm('create', cfg_Recon);
          
          op.sino_id = sino_id;
          op.rec_id  = rec_id;
          op.proj_size = proj_size;
          op.rec_size  = vol_size;
          op.rec_id = rec_id;
%           op.model_type = type;
          op.proj_geo   = proj_geo;
%           op.rec_handle  = rec_handle;
          op.cflag       = false;
          op.sweepflag   = false;
          op.num_iter = numiter ;
          op.funHandle = @opRecon_multiply;
           
       end
   % delete functions are necessary to properly delete ASTRA structures from memory, otherwise there will be a
   % memory leak and after creating several operators matlab will crash and give "out of memory" error
   % Delete may have only 1 argument, and since for different structures different type of 'ASTRA delete' is
   % necessary there should be 3 delete functions: for sinogram, for astra_data_2d and Recon_operator 
%        function delete(op)
%            %...
%        end

   end
   
   methods( Access = protected )

        function y = multiply(op,x,mode)
            y = op.funHandle(op, x, mode);
        end

   end
   
   
   methods( Access = private )

        % 2D projection code
        function y = opRecon_multiply(op,x,mode)
            
            if mode == 1
                if ~isempty(op.x0)
                    if size(op.x0,2) == 1
                        op.x0 = reshape(op.x0, op.rec_size);
                    end
                    astra_mex_data2d('store', op.rec_id, op.x0);
                end
                x = reshape(x, op.proj_size);
                astra_mex_data2d('store', op.sino_id, x);
                
                astra_mex_algorithm('iterate', op.alg_id, op.num_iter);
                y = astra_mex_data2d('get', op.rec_id);              %%%%%%!!!!!!!!!!!!! mogê siê tak odwo³aæ?
                y = y(:);

            else
                error('Inverse operator not implemented')
            end
        end
   end
end