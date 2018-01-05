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
params_ini=[find(mean_1 == max(mean_1)) max(mean_1) find(mean_1>=0.4,1) 0.5
    find(mean_2 == max(mean_2)) max(mean_2) find(mean_2>=0.4,1) 0.5
    find(mean_3 == max(mean_3)) max(mean_3) find(mean_3>=0.4,1) 0.5];

a=zeros(N,4);
gammatheo=zeros(size(i_pixels));

file = fopen(['params_' ind_patient '.txt'],'w');
fprintf(file , [  'tmax' '\t' 'ymax' '\t' 'd' '\t' 'a' '\r\n']);
for i=1:N
    disp(i)
    if max(i_pixels(i,:))==0
        a(i,:)=[-1 0 -1 -1];
        gammatheo(i,:)=zeros(1,tmax);
    else
        a(i,:)=levenberg_marquardt(1:tmax, i_pixels(i,:), params_ini(m_pixels(i,4),:), 0.001, 10000);
        gammatheo(i,:)=arrayfun(@(x) gamma(x,a(i,:)),1:tmax);
    end
    fprintf(file,'%f\t',a(i,:));
    fprintf(file, '\r\n');
end
fclose(file);

R2 = sum(i_pixels-mean(i_pixels,2),2)./sum(i_pixels-gammatheo,2);

% a=levenberg_marquardt(1:tmax, i_pixels(i,:), [find(mean_1 == max(mean_1)) max(mean_1) find(mean_1>=0.4,1) 0.5], 0.001, 10000)
% gammatheo=arrayfun(@(x) gamma(x,a),1:tmax);
% 
% nbfig=length(findobj('type','figure'));
% figure(nbfig+1);
% plot(1:tmax, mean_1, 'r', 1:tmax, mean_2, 'b', 1:tmax, mean_3, 'g',1:tmax,gammatheo,'m')
% xlabel('Temps')
% ylabel('Concentration')
% legend('Tissus Gris','Tissus Blancs','L�sion','Gamma theo (gris)')

end


function f=gamma(t, a)
    tmax = a(1);
    ymax = a(2);
    d = a(3);
    alpha = a(4);
    if (t<=d)
        f = 0;
    else
        T = (t-d)/tmax;
        f = ymax*T^alpha*exp(alpha*(1-T));
    end
end