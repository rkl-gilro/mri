%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the folder where the results are
folder = '/home/raquel/Documents/repositories/mri-segmentation/resources/registration/resultsT1/lassalas/nonrigiddeformation/deformation3/';
dcm = dir(folder);

%% Open a txt file
fileName = strcat(folder,'/results_output.txt');
[fid, message] = fopen(fileName,'a');


for i=3:size(dcm,1)-1
        
    filename = strcat(folder,dcm(i).name,'/lassalas.mat');
    
    %% Load the matlab file
    load(filename);
    
    if deformation_type.rigid 
        
        trans1 = sum(abs(deformation_type.translation_ax));
        trans2 = sum(abs(deformation_type.translation_sag));
        trans3 = sum(abs(deformation_type.translation_cor));
        
        ang1 = sum(abs(deformation_type.angles_ax));
        ang2 = sum(abs(deformation_type.angles_sag));
        ang3 = sum(abs(deformation_type.angles_cor));
        
        %% Write results in txt file
        fprintf(fid,'%d \t %5.3f \t %5.3f \t %5.3f \t %5.3f \t %5.3f \t %5.3f \t %5.3f + %5.3f \t %5.3f + %5.3f \t %5.3f + %5.3f \n', ... %5.3f, %5.3f \t %5.3f , %5.3f \t %5.3f , %5.3f \n', ...
                 i-2, trans1, trans2, trans3, ang1, ang2, ang3, mean(error_int1(:)), std(error_int1(:)),mean(error_int2(:)), std(error_int2(:)),mean(error_int3(:)), std(error_int3(:)));
             
%                  mean(straight_ax_int1(:)), mean(straight_sag_int1(:)), mean(straight_ax_int2(:)), mean(straight_cor_int2(:)), mean(straight_cor_int3(:)), mean(straight_sag_int3(:)));
             
    end
    
    if deformation_type.nonrigid
        
        magnitude = sum(abs(deformation_type.current_mag_def));
        
        %% Write results in txt file
        fprintf(fid,'%d \t %5.3f \t %5.3f + %5.3f \t %5.3f + %5.3f \t %5.3f + %5.3f \n', ...
                 i-2, magnitude, mean(error_int1(:)), std(error_int1(:)),mean(error_int2(:)), std(error_int2(:)),mean(error_int3(:)), std(error_int3(:)));
%                  mean(straight_ax_int1(:)), mean(straight_sag_int1(:)), mean(straight_ax_int2(:)), mean(straight_cor_int2(:)), mean(straight_cor_int3(:)), mean(straight_sag_int3(:)));
             
    end
    
    

end

fclose(fid);

