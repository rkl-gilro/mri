%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main reading patient from a folder with DICOMDIR
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

start_path = 'resources/patients/';
folder_name = uigetdir(start_path);

%% Read the dicomdir and save the info in dicomdir.mat in 'folder'
load_dicomdir(folder_name);

mat_file = strcat(folder_name,'/dicomdir.mat');
load(mat_file);

for i=1:length(dcmdir.dcmPatient.Study.Series)
    series_des{i}.images = dcmdir.dcmPatient.Study.Series(i,1).Images;
    series_des{i}.info   = dcmdir.dcmPatient.Study.Series(i,1).ImagesInfo;
    name_series{i} = dcmdir.dcmPatient.Study.Series(i,1).SeriesDescription;
end
%% Define a struct with axial, sagittal and coronal views
% views = get_3views(dcmdir);
views = choose_ax_sag_cor(series_des, name_series);

%% Clear the rest of the data

% clearvars -except views dcmdir
% clc

