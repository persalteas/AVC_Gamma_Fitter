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

%Matrice des masques
masque048 = load_untouch_nii(input_mask048);
masque249 = load_untouch_nii(input_mask249);
masque260 = load_untouch_nii(input_mask260);
masque170 = load_untouch_nii(input_mask170);
% calcul des pixels parenchymateux
masque170.img(:,:,17);

