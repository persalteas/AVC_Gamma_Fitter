% ind_patient :  (char)   index du patient

function AnalysePatient(ind_patient)
disp(['Analyse patient ' ind_patient])
%%%%%%%%%%%%%%%%% OUVERTURE DES MASQUES %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_folder = [ pwd '\Masques\'];

% Lecture du masque
input_file = [ mask_folder 'mask_patient_' ind_patient '.nii' ];
mask = load_untouch_nii(input_file);

% Matrice du masque
matrix_m = mask.img(:,:,:);
%matrix_m = matrix_m(:,:,1:20); % A supprimer si utilisation des bons masques

disp('Extraction des pixels parenchymateux des masques...')

%Creation fichier txt pour pixels parenchymateux 
file = fopen(['mask_' ind_patient '.txt'],'w');
fprintf(file , [  'Dim1' '\t' 'Dim2' '\t' 'Coupe' '\t' 'Value' '\r\n']);

%Creer un vecteur des pixels parenchymateux
N = sum(sum(sum(matrix_m ~= 0))); % Nombre de pixels parenchymateux
m_pixels = zeros(N, 4);

i = 1;
for coupe=1:1:size(matrix_m,3)
    for dim1=1:1:size(matrix_m,1)
        for dim2=1:1:size(matrix_m,2)
            if(mask.img(dim1,dim2,coupe)>=1)
               fprintf(file, '%d\t%d\t%d\t%f\r\n', [dim1 dim2 coupe mask.img(dim1,dim2,coupe)]');
               m_pixels(i,:) = [ dim1 dim2 coupe mask.img(dim1,dim2,coupe)];
               i = i+1;
            end 
        end
    end
end
fclose(file);  

%%%%%%%%%%%%%%%%%%% OUVERTURE DES IMAGES %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_folder = [ pwd '\Images_Concentration\'];

% Lecture de l'image
input_file = [ image_folder 'ConcentrationImage_patient_048.nii' ];
pwi = load_untouch_nii(input_file);

% Matrice de l'image
matrix_i = pwi.img(:,:,:,:);

disp('Extraction des pixels parenchymateux des images...')

% Creer une matrice des pixels parenchymateux
tmax = size(matrix_i,4); 
i_pixels = zeros(N, tmax);  

% Creation fichier txt pour pixels parenchymateux 

file = fopen(['img_' ind_patient '.txt'],'w');
for i=1:N
    x = m_pixels(i,1);
    y = m_pixels(i,2);
    z = m_pixels(i,3);
    i_pixels(i,:) = reshape(matrix_i(x,y,z,1:tmax),1,tmax);
    fprintf(file,'%f\t',i_pixels(i,:));
    fprintf(file, '\r\n');
end
i_pixels(i_pixels<0) = 0;
fclose(file);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Affichage des dynamiques 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_1 = mean(i_pixels(m_pixels(:,4)==1,:));
mean_2 = mean(i_pixels(m_pixels(:,4)==2,:));
mean_3 = mean(i_pixels(m_pixels(:,4)==3,:));

a=levenberg_marquardt(1:tmax, mean_1, [25.0 6.0 10.0 1.0], 0.001, 10000);
a

nbfig=length(findobj('type','figure'));
figure(nbfig+1);
plot(1:tmax, mean_1, 'r', 1:tmax, mean_2, 'b', 1:tmax, mean_3, 'g')
xlabel('Temps')
ylabel('Concentration')
legend('Tissus Gris','Tissus Blancs','Lésion')

end