function views = get_3views(dicom_series)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  From the given struct 'dicom_series', it extracts 
%%  the 3 views (axial, sagittal and coronal)
%% 
%%  Inputs:  1. dicom_series -> struct containing dicomdir info
%%  Outputs: 1. views -> struct containing axial, sagittal, coronal
%%                       matrices of size (512,512,s), and its corresponding
%%                       dicom_info
%%                      e.g., views.axial & views.axial_info
%%
%% Execute:
%% -folder = 'resources/patients/patient1/' % folder containing dicomdir file
%% -load_dicomdir(folder);
%% -views = get_3views(dcmdir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Most of the times the second series corresponds to sagittal

for i = 1:length(dicom_series.dcmPatient.Study.Series)
    
    if strcmp(dicom_series.dcmPatient.Study.Series(i,1).SeriesDescription,'SAG T2 PROPELLER')
        series_sag = i;
    elseif strcmp(dicom_series.dcmPatient.Study.Series(i,1).SeriesDescription,'Ax T2 PROPELLER')
        series_ax = i;
    elseif strcmp(dicom_series.dcmPatient.Study.Series(i,1).SeriesDescription,'COR T2 PROPELLER')
        series_cor = i;
    end
    
end


%% Series Sagittal
[m n] = size(dicom_series.dcmPatient.Study.Series(series_sag,1).Images{1});
[p q] = size(dicom_series.dcmPatient.Study.Series(series_sag,1).Images);

for i = 1:p
    views.sagittal(:,:,i) = dicom_series.dcmPatient.Study.Series(series_sag,1).Images{i};
    views.sagittal_info{i} = dicom_series.dcmPatient.Study.Series(series_sag,1).ImagesInfo{i};
end
%% Series Coronal
[m n] = size(dicom_series.dcmPatient.Study.Series(series_cor,1).Images{1});
[p q] = size(dicom_series.dcmPatient.Study.Series(series_cor,1).Images);

for i = 1:p
    views.coronal(:,:,i) = dicom_series.dcmPatient.Study.Series(series_cor,1).Images{i};
    views.coronal_info{i} = dicom_series.dcmPatient.Study.Series(series_cor,1).ImagesInfo{i};
end

%% Series Axial
[m n] = size(dicom_series.dcmPatient.Study.Series(series_ax,1).Images{1});
[p q] = size(dicom_series.dcmPatient.Study.Series(series_ax,1).Images);

for i = 1:p
    views.axial(:,:,i) = dicom_series.dcmPatient.Study.Series(series_ax,1).Images{i};
    views.axial_info{i} = dicom_series.dcmPatient.Study.Series(series_ax,1).ImagesInfo{i};
end

% dcmdir.dcmPatient.Study.Series(i,1).SeriesDescription
