image_folder = 'Images_Concentration\';
mask_folder = 'Masques\';

input_file = [ image_folder 'ConcentrationImage_patient_048.nii' ];
pwi_048 = load_untouch_nii(input_file);

input_file = [ image_folder 'ConcentrationImage_patient_249.nii' ];
pwi_249 = load_untouch_nii(input_file);

input_file = [ image_folder 'ConcentrationImage_patient_260.nii' ];
pwi_260 = load_untouch_nii(input_file);

input_file = [ image_folder 'ConcentrationImage_patient_170.nii' ];
pwi_170 = load_untouch_nii(input_file);


input_file = [ mask_folder 'mask_patient_048.nii' ];
mask_048 = load_untouch_nii(input_file);

input_file = [ mask_folder 'mask_patient_249.nii' ];
mask_249 = load_untouch_nii(input_file);

input_file = [ mask_folder 'mask_patient_260.nii' ];
mask_260 = load_untouch_nii(input_file);

input_file = [ mask_folder 'mask_patient_170.nii' ];
mask_170 = load_untouch_nii(input_file);

%Matrice des masques complètes
matrix_048 = mask_048.img(:,:,:);
matrix_249 = mask_249.img(:,:,:);
matrix_260 = mask_260.img(:,:,:);
matrix_170 = mask_170.img(:,:,:);

%Dim matrice, Nombre de coupes
size(matrix_048)
size(matrix_249)
size(matrix_260)
size(matrix_170)


%Creation fichier text pour pixels parenchymateux pour les quatres masques
file048 = fopen('mask_048.txt','w');
fprintf(file048 , [ 'Coupe' '\t' 'Dim1' '\t' 'Dim2' '\t' 'Value' '\r\n']);

file249 = fopen('mask_249.txt','w');
fprintf(file249 , [ 'Coupe' '\t' 'Dim1' '\t' 'Dim2' '\t' 'Value' '\r\n']);

file260 = fopen('mask_260.txt','w');
fprintf(file260 , [ 'Coupe' '\t' 'Dim1' '\t' 'Dim2' '\t' 'Value' '\r\n']);

file170 = fopen('mask_170.txt','w');
fprintf(file170 , [ 'Coupe' '\t' 'Dim1' '\t' 'Dim2' '\t' 'Value' '\r\n']);

for coupe=1:1:size(matrix_048,3)
    for dim1=1:1:size(matrix_048,1)
        for dim2=1:1:size(matrix_048,2)
            if(mask_048.img(dim1,dim2,coupe)~=0)
               fprintf(file048, '%d \t %d \t %d %f \r\n', [coupe dim1 dim2 mask_048.img(dim1,dim2,coupe)]');
            end 
            if(mask_249.img(dim1,dim2,coupe)~=0)
               fprintf(file249, '%d \t %d \t %d %f \r\n', [coupe dim1 dim2 mask_249.img(dim1,dim2,coupe)]');
            end 
            if(mask_260.img(dim1,dim2,coupe)~=0)
               fprintf(file260, '%d \t %d \t %d %f \r\n', [coupe dim1 dim2 mask_260.img(dim1,dim2,coupe)]');
            end 
            if(mask_170.img(dim1,dim2,coupe)~=0)
               fprintf(file170, '%d \t %d \t %d %f \r\n', [coupe dim1 dim2 mask_170.img(dim1,dim2,coupe)]');
            end 
        end
    end
end
fclose(file048);  
fclose(file249);  
fclose(file260);  
fclose(file170);  



