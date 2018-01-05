clear all;
addpath('NIfTI_20140122');

image_folder = [ pwd '\Images_Concentration\'];
mask_folder = [ pwd '\Masques\'];

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
matrix_m_048 = mask_048.img(:,:,:);
matrix_m_249 = mask_249.img(:,:,:);
matrix_m_260 = mask_260.img(:,:,:);
matrix_m_170 = mask_170.img(:,:,:);
matrix_i_048 = pwi_048.img(:,:,:,:);
matrix_i_249 = pwi_249.img(:,:,:,:);
matrix_i_260 = pwi_260.img(:,:,:,:);
matrix_i_170 = pwi_170.img(:,:,:,:);

% matrix_m_048 = matrix_m_048(:,:,1:20);
% matrix_m_249 = matrix_m_249(:,:,1:20);
% matrix_m_260 = matrix_m_260(:,:,1:20);
% matrix_m_170 = matrix_m_170(:,:,1:20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUVERTURE DES MASQUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extraction des pixels parenchymateux des masques...')


%Dim matrice, Nombre de coupes
size(matrix_m_048);
size(matrix_m_249);
size(matrix_m_260);
size(matrix_m_170);

%Creation fichier text pour pixels parenchymateux pour les quatres masques
file048 = fopen('mask_048.txt','w');
fprintf(file048 , [  'Dim1' '\t' 'Dim2' '\t' 'Coupe' '\t' 'Value' '\r\n']);

file249 = fopen('mask_249.txt','w');
fprintf(file249 , [  'Dim1' '\t' 'Dim2' '\t' 'Coupe' '\t' 'Value' '\r\n']);

file260 = fopen('mask_260.txt','w');
fprintf(file260 , [  'Dim1' '\t' 'Dim2' '\t' 'Coupe' '\t' 'Value' '\r\n']);

file170 = fopen('mask_170.txt','w');
fprintf(file170 , [  'Dim1' '\t' 'Dim2' '\t' 'Coupe' '\t' 'Value' '\r\n']);

%Creer un vecteur des pixels parenchymateux
N48 = sum(sum(sum(matrix_m_048 ~= 0)));
N249 = sum(sum(sum(matrix_m_249 ~= 0)));
N260 = sum(sum(sum(matrix_m_260 ~= 0)));
N170 = sum(sum(sum(matrix_m_170 ~= 0)));

m_pixels48 = zeros(N48, 4);
m_pixels249 = zeros(N249, 4);
m_pixels260 = zeros(N260, 4);
m_pixels170 = zeros(N170, 4);

i = 1;
j = 1;
k = 1;
l = 1;
for coupe=1:1:size(matrix_m_048,3)
    for dim1=1:1:size(matrix_m_048,1)
        for dim2=1:1:size(matrix_m_048,2)
            if(mask_048.img(dim1,dim2,coupe)>=1)
               fprintf(file048, '%d\t%d\t%d\t%f\r\n', [dim1 dim2 coupe mask_048.img(dim1,dim2,coupe)]');
               m_pixels48(i,:) = [ dim1 dim2 coupe mask_048.img(dim1,dim2,coupe)];
               i = i+1;
            end 
            if(mask_249.img(dim1,dim2,coupe)>=1)
               fprintf(file249, '%d\t%d\t%d\t%f\r\n', [ dim1 dim2 coupe mask_249.img(dim1,dim2,coupe)]');
               m_pixels249(j,:) = [ dim1 dim2 coupe mask_249.img(dim1,dim2,coupe)];
               j = j+1;
            end
            if(mask_260.img(dim1,dim2,coupe)>=1)
               fprintf(file260, '%d\t%d\t%d\t%f\r\n', [ dim1 dim2 coupe mask_260.img(dim1,dim2,coupe)]');
               m_pixels260(k,:) = [ dim1 dim2 coupe mask_260.img(dim1,dim2,coupe)];
               k = k+1;
            end 
            if(mask_170.img(dim1,dim2,coupe)>=1)
               fprintf(file170, '%d\t%d\t%d\t%f\r\n', [ dim1 dim2 coupe mask_170.img(dim1,dim2,coupe)]');
               m_pixels170(l,:) = [ dim1 dim2 coupe mask_170.img(dim1,dim2,coupe)];
               l = l+1;
            end 
        end
    end
end
fclose(file048);  
fclose(file249);  
fclose(file260);  
fclose(file170);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIXELS DES IMAGES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extraction des pixels parenchymateux des images...')
tmax = size(matrix_i_048,4);  % c a d 60
i_pixels48 = zeros(N48, tmax);  
i_pixels249 = zeros(N249, tmax);
i_pixels260 = zeros(N260, tmax);
i_pixels170 = zeros(N170, tmax);


% %Creation fichier text pour pixels parenchymateux pour les quatres masques

file048 = fopen('img_048.txt','w');
for i=1:N48
    x = m_pixels48(i,1);
    y = m_pixels48(i,2);
    z = m_pixels48(i,3);
    i_pixels48(i,:) = reshape(matrix_i_048(x,y,z,1:tmax),1,tmax);
    fprintf(file048,'%f\t',i_pixels48(i,:));
    fprintf(file048, '\r\n');
end
i_pixels48(i_pixels48<0) = 0;
fclose(file048);

file249 = fopen('img_249.txt','w');
for i=1:N249
    x = m_pixels249(i,1);
    y = m_pixels249(i,2);
    z = m_pixels249(i,3);
    i_pixels249(i,:) = reshape(matrix_i_249(x,y,z,1:tmax),1,tmax);
    fprintf(file249,'%f\t',i_pixels249(i,:));
    fprintf(file249, '\r\n');
end
fclose(file249);

file260 = fopen('img_260.txt','w');
for i=1:N260
    x = m_pixels260(i,1);
    y = m_pixels260(i,2);
    z = m_pixels260(i,3);
    i_pixels260(i,:) = reshape(matrix_i_260(x,y,z,1:tmax),1,tmax);
    fprintf(file260,'%f\t',i_pixels260(i,:));
    fprintf(file260, '\r\n');
end
fclose(file260);

file170 = fopen('img_170.txt','w');
for i=1:N170
    x = m_pixels170(i,1);
    y = m_pixels170(i,2);
    z = m_pixels170(i,3);
    i_pixels170(i,:) = reshape(matrix_i_170(x,y,z,1:tmax),1,tmax);
    fprintf(file170,'%f\t',i_pixels170(i,:));
    fprintf(file170, '\r\n');
end
fclose(file170);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Affichage des dynamiques 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean48_1 = mean(i_pixels48(m_pixels48(:,4)==1,:));
mean48_2 = mean(i_pixels48(m_pixels48(:,4)==2,:));
mean48_3 = mean(i_pixels48(m_pixels48(:,4)==3,:));

hold on;
figure(1);
plot(1:tmax, mean48_1, 'r', 1:tmax, mean48_2, 'b', 1:tmax, mean48_3, 'g')
hold off;