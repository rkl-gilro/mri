%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Main read a folder that contains DICOM files
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

start_path = 'resources/patients/';
folder_name = uigetdir(start_path);

mat_file = strcat(folder_name,'/dcm.mat');

if exist(mat_file,'file')
    tmp = load(mat_file);
    series_des  = tmp.series_des;
    name_series = tmp.name_series;
else
    [~, ~, series_des, name_series] = read_dicom(folder_name, 1, 0);
    save(mat_file, 'series_des', 'name_series');
end

views = choose_ax_sag_cor(series_des, name_series);

clearvars -except views series_des name_series