function create_dicom( images, info, save_name )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%% Inputs:  1. images    -> (cell) or (matrix) containing the images
%%          2. info      -> (cell) that contains the dicom_info of the images
%%          3. save_name -> (string) name of the folder to save the dicom
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 0;

%% Check the input arguments
if iscell(images)
    
    for k=1:length(images)

        for i=1:size(images{k},3)
            
            a = i + b;
            
            [~, add] = add_zeros(a, 3);
            
            file_name = strcat(save_name,'mri_file',add, num2str(a),'.dcm');
            dicomwrite(int16(images{k}(:,:,i)), file_name, info{k}{i});
        end
        b = b + size(images{k},3);
    end
    
else
    
    for i=1:size(images,3)
        file_name = strcat(save_name,'mri_file',num2str(i),'.dcm');
        dicomwrite(int16(images(:,:,i)), file_name, info{i});
    end
    
end