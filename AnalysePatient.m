% ind_patient :  (char)   index du patient

function AnalysePatient(ind_patient, display)
	disp(['Analyse patient ' ind_patient])

	%%%%%%%%%%%%%%%%% OUVERTURE DES MASQUES %%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	mask_folder = [ pwd '\Masques\'];

	% Lecture du masque
	input_file = [ mask_folder 'mask_patient_' ind_patient '.nii' ];
	mask = load_untouch_nii(input_file);

	% Matrice du masque
	matrix_m = mask.img(:,:,:);

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
					fprintf(file, '%d\t%d\t%d\t%f\r\n', [dim1 dim2 coupe mask.img(dim1,dim2,coupe)]);
					m_pixels(i,:) = [ dim1 dim2 coupe mask.img(dim1,dim2,coupe) ];
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
	tmax = tmax - 10; % Suppression arbitraire des 10 derniers (artefact qui ne colle pas à la loi gamma)
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

	% %%%%%%%%%%%%%%%% Estimation des parametres %%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Analyse des pixels...')

	mean_1 = mean(i_pixels(m_pixels(:,4)==1,:)); % Moyenne d'intensite des pixels parenchymateux de classe 1
	mean_2 = mean(i_pixels(m_pixels(:,4)==2,:)); % Moyenne d'intensite des pixels parenchymateux de classe 2
	mean_3 = mean(i_pixels(m_pixels(:,4)==3,:)); % Moyenne d'intensite des pixels parenchymateux de classe 3

	% Valeurs initiales des parametres a optimiser, selon la classe de pixel
	% tmax = abscisse du pic
	% ymax = hauteur du pic
	% d = ? (pourquoi 0.4 ?)
	% alpha = parametre de forme de la fonction gamma
	%			[	  tmax						ymax			d		   alpha ]
	params_ini = [  find(mean_1 == max(mean_1)) max(mean_1) find(mean_1>=0.4,1) 2.0
					find(mean_2 == max(mean_2)) max(mean_2) find(mean_2>=0.4,1) 2.0
					find(mean_3 == max(mean_3)) max(mean_3) find(mean_3>=0.4,1) 2.0   ];

	params = zeros(N,4);
	gammatheo = zeros(size(i_pixels));

	file = fopen(['params_' ind_patient '.txt'],'w');
	fprintf(file , [  'tmax' '\t' 'ymax' '\t' 'd' '\t' 'a' '\r\n']);
	for i=1:N
		if (display==1 || display==3), fprintf('analyse du pixel %d... ',i), end
		if max(i_pixels(i,:)) == 0 % ce pixel est toujours eteint
			params(i,:) = [-1 0 -1 -1];		  % params
			gammatheo(i,:) = zeros(1,tmax); % shape of the gamma distribution
			if (display==1 || display==3), disp('pixel toujours eteint.'), end
		else
			params(i,:) = levenberg_marquardt(1:tmax, i_pixels(i,:), params_ini(m_pixels(i,4),:), 0.1, 10000, i, display); % params
			gammatheo(i,:) = arrayfun(@(x) gamma(x,params(i,:)),1:tmax); % shape of the gamma distribution
		end
		fprintf(file,'%f\t',params(i,:));
		fprintf(file, '\r\n');
	end
	fclose(file);
	fprintf('Parametres des lois Gamma estimees sauvegardes dans params_%s.txt\n', ind_patient);
	R2 = sum(i_pixels-mean(i_pixels,2),2)./sum(i_pixels-gammatheo,2);
    R2_mat=ones(size(matrix_m))*NaN;
    M=2;
    m=0;
    for i = 1:N
        if(R2(i)<0)
           R2_mat(m_pixels(i,1),m_pixels(i,2), m_pixels(i,3))= 0;
        else if (R2(i)>2)
           R2_mat(m_pixels(i,1),m_pixels(i,2), m_pixels(i,3))= 255;
            else       
                R2_mat(m_pixels(i,1),m_pixels(i,2), m_pixels(i,3))=R2(i)*255/(M-m)+255*m/(m-M);
            end
        end 
    end
end

% la fonction de distribution gamma (fonction du temps)
function f=gamma(t, params)
	tmax = params(1);
	ymax = params(2);
	d = params(3);
	alpha = params(4);
	if (t<=d)
		f = 0;
	else
		T = (t-d)/tmax;
		f = ymax * T^alpha * exp(alpha*(1-T));
	end
end