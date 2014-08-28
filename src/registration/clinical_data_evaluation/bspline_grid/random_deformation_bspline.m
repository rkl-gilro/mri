function [vol_def, grid_def] = random_deformation_bspline( vol, max_def, percent )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Apply a non rigid deformation using b-spline grid 
%% 'D.Kroon University of Twente (February 2010)'
%% b-spline-grid--image-and-point-based-registration
%%
%% Inputs:  1. vol -> MxNxS image volume
%%          2. max_def -> scalar that gives the magnitude of deformation
%%          3. percent -> scalar [0, 1], 0-take the hold volume for the deformation
%%                                       1-no deformation is applied 
%%             the smaller the value, the bigger portion of volume is taken for the deformation
%%
%% Outputs: 1. vol_reg -> MxNxS deformed image volume
%%          2. grid_def -> struct containing 
%%                         2.1. size -> 3x1 array with the size of the grid 
%%                         2.2. xx, yy, zz -> array the coordinates of the original grid 
%%                         2.3. xx_def, yy_def, zz_def -> array of coordinates of deformed grid
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


defaultoptions = struct('Similarity',[],'Registration','Both','Penalty',1e-3,'MaxRef',2,'Grid',[],'Spacing',[1 1 1],'MaskMoving',[],'MaskStatic',[],'Verbose',2,'Points1',[],'Points2',[],'PStrength',[],'Interpolation','Linear','Scaling',[1 1]);
Options = defaultoptions;

M = eye(4);

% Make the initial b-spline registration grid (uniform)
[O_trans, Spacing, ~] = Make_Initial_Grid(Options.Grid, Options.Spacing, Options, vol, M, true);

xx = O_trans(:,:,:,1);
xx = xx(:);
yy = O_trans(:,:,:,2);
yy = yy(:);
zz = O_trans(:,:,:,3);
zz = zz(:);

% Compute the random displacements for each direction
a = -max_def; b = max_def; 

sigma = b-a;

nx = size(O_trans(:,:,:,1),1);
ny = size(O_trans(:,:,:,1),2);
nz = size(O_trans(:,:,:,1),3);

out_ind = ones(nx * ny * nz, 1);

% Define the percentage of volume that is gonna be used for applying the
% deformation

range_x = percent*(max(xx) - min(xx))/2;
range_y = percent*(max(yy) - min(yy))/2;
range_z = percent*(max(zz) - min(zz))/2;

% Get the indexes of the grid where the deformation is gonna be applied
for i = 1:nx
    
    for j = 1:ny
        
        for tmp =  1:nz
        
            s2ind =  tmp + nz*(j-1 + ny*(i-1));
            

            if xx(s2ind) >= min(xx) + range_x && xx(s2ind) <= max(xx) - range_x && yy(s2ind) >= min(yy) + range_y && yy(s2ind) <= max(yy) - range_y && zz(s2ind) >= min(zz) + range_z && zz(s2ind) <= max(zz) - range_z % O_trans(i,j,tmp,1) >= 100 && O_trans(i,j,tmp,1) <= 400 && O_trans(i,j,tmp,2) >= 100 && O_trans(i,j,tmp,2) <= 400 %(i>=thr && i<=nx-thr) && (j>=thr && j<= ny-thr) && (tmp>=1 && tmp<=nz-1)  %  i==1 || i==nx || j==1 || j==ny || tmp==1 || tmp==nz 
                    
               out_ind(s2ind) = 0;
                
            end

        end
    end
end

% Apply random displacements to the grid
xx_def(out_ind==1) = xx(out_ind==1); 
yy_def(out_ind==1) = yy(out_ind==1); 
zz_def(out_ind==1) = zz(out_ind==1);

xx_def(out_ind==0) = xx(out_ind==0) + a + sigma.*rand(size(xx(out_ind==0)));
yy_def(out_ind==0) = yy(out_ind==0) + a + sigma.*rand(size(yy(out_ind==0)));
zz_def(out_ind==0) = zz(out_ind==0) -.5 + 1.*rand(size(zz(out_ind==0)));

% Define the grid again
O_trans(:,:,:,1) = reshape(xx_def, size(O_trans(:,:,:,1)));
O_trans(:,:,:,2) = reshape(yy_def, size(O_trans(:,:,:,2)));
O_trans(:,:,:,3) = reshape(zz_def, size(O_trans(:,:,:,3)));

% Define the output
grid_def.size = size(O_trans(:,:,:,1));
grid_def.spacing = Spacing;

grid_def.x = xx;
grid_def.y = yy;
grid_def.z = zz;

grid_def.xx_def = xx_def;
grid_def.yy_def = yy_def;
grid_def.zz_def = zz_def;

% Calculate the new volume
vol_def = bspline_transform(O_trans, vol, Spacing, 3);


function [O_trans, Spacing, MaxItt] = Make_Initial_Grid(O_trans, Spacing, Options, Imoving, M, IS3D)
if(isempty(O_trans)),
    
    if(isempty(Options.Spacing))
        if(IS3D)
            % Calculate max refinements steps
            MaxItt=min(floor(log2(size(Imoving)/4)));
            
            % set b-spline grid spacing in x,y and z direction
            Spacing=[2^MaxItt 2^MaxItt 2^MaxItt];
        else
            % Calculate max refinements steps
            MaxItt=min(floor(log2([size(Imoving,1) size(Imoving,2)]/4)));
            
            % set b-spline grid spacing in x and y direction
            Spacing=[2^MaxItt 2^MaxItt];
        end
    else
        % set b-spline grid spacing in x and y direction
        Spacing=round(Options.Spacing);
        t=Spacing; MaxItt=0; while((nnz(mod(t,2))==0)&&(nnz(t<8)==0)), MaxItt=MaxItt+1; t=t/2; end
    end
    
    Spacing = [20 20 4];
    
    % Make the Initial b-spline registration grid
    if(IS3D)
        O_trans=make_init_grid(Spacing,size(Imoving),M);
    else
        if(strcmpi(Options.Registration(1),'N'))
            O_trans=make_init_grid(Spacing,[size(Imoving,1) size(Imoving,2)]);
        else
            O_trans=make_init_grid(Spacing,[size(Imoving,1) size(Imoving,2)],M);
        end
    end
else
    MaxItt=0;
    TestSpacing=Spacing;
    while(mod(TestSpacing,2)==0), TestSpacing=TestSpacing/2; MaxItt=MaxItt+1; end
    
    if(IS3D)
        % Calculate center of the image
        mean=size(Imoving)/2;
        % Make center of the image coordinates 0,0
        xd=O_trans(:,:,:,1)-mean(1); yd=O_trans(:,:,:,2)-mean(2); zd=O_trans(:,:,:,3)-mean(3);
        % Calculate the rigid transformed coordinates
        O_trans(:,:,:,1) = mean(1) + M(1,1) * xd + M(1,2) *yd + M(1,3) *zd + M(1,4)* 1;
        O_trans(:,:,:,2) = mean(2) + M(2,1) * xd + M(2,2) *yd + M(2,3) *zd + M(2,4)* 1;
        O_trans(:,:,:,3) = mean(3) + M(3,1) * xd + M(3,2) *yd + M(3,3) *zd + M(3,4)* 1;
    else
        if(~strcmpi(Options.Registration(1),'N'))
            % Calculate center of the image
            mean=size(Imoving)/2;
            % Make center of the image coordinates 0,0
            xd=O_trans(:,:,1)-mean(1); yd=O_trans(:,:,2)-mean(2);
            % Calculate the affine transformed coordinates
            O_trans(:,:,1) = mean(1) + M(1,1) * xd + M(1,2) *yd + M(1,3) * 1;
            O_trans(:,:,2) = mean(2) + M(2,1) * xd + M(2,2) *yd + M(2,3) * 1;
        end
    end
end
% Limit refinements steps to user input
if(Options.MaxRef<MaxItt), MaxItt=Options.MaxRef; end

