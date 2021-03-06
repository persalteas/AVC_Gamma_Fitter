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
	tmax = tmax - 10; % Suppression arbitraire des 10 derniers (artefact qui ne colle pas � la loi gamma)
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
			params(i,:) = [NaN NaN NaN NaN];		  % params
			%gammatheo(i,:) = zeros(1,tmax); % shape of the gamma distribution
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

    % Calcul du R2
    Max_classes=zeros(N,1);
    Max_classes(m_pixels(:,4)==1)=max(max(i_pixels(m_pixels(:,4)==1,:),[],2));
    Max_classes(m_pixels(:,4)==2)=max(max(i_pixels(m_pixels(:,4)==2,:),[],2));
    Max_classes(m_pixels(:,4)==3)=max(max(i_pixels(m_pixels(:,4)==3,:),[],2));
    R2 = sum(abs(i_pixels-gammatheo),2)./(Max_classes.*tmax);

    % Convertion en matrice pour l'afficher
    R2_mat=zeros(size(matrix_m));
    M=max(R2);
    m=min(R2);
    for i = 1:N
        R2_mat(m_pixels(i,1),m_pixels(i,2), m_pixels(i,3))=R2(i)*255/(M-m)+255*m/(m-M);
    end
    
    % Affichage du R2 par coupe
     for coupe = 1:size(R2_mat,3)
         image(R2_mat(:,:,coupe),'CDataMapping','scaled')
         title(['coupe : ' coupe])
         colorbar
         hold off;
         waitforbuttonpress; % Attend le clic de l'utilisateur sur le graphique pour passer au suivant
         close;
     end    
    
    % Sauvegarde au format nii
    nii=make_nii(R2_mat);
    save_nii(nii,['result_' ind_patient '.nii'])
    
    % Statistiques sur les param�tres
    mean(params(m_pixels(:,4)==1,:),'omitnan')
    mean(params(m_pixels(:,4)==2,:),'omitnan')
    mean(params(m_pixels(:,4)==3,:),'omitnan')
    std(params(m_pixels(:,4)==1,:),'omitnan')
    std(params(m_pixels(:,4)==2,:),'omitnan')
    std(params(m_pixels(:,4)==3,:),'omitnan')
    
    %Tmax
    figure(length(findobj('type','figure'))+1) 
    title(['Patient ' num2str(ind_patient)])
    for i=1:3
        subplot(3,1,i)
        p_tmax=params(m_pixels(:,4)==i,1);
        %hist(p_tmax(p_tmax<50))
        %title(['tmax ' num2str(ind_patient)])
        histogram(p_tmax(p_tmax<50),15)
        xlim([0,50])
    end
    
    %Ymax
    figure(length(findobj('type','figure'))+2)
    for i=1:3
        subplot(3,1,i)
        p_ymax=params(m_pixels(:,4)==i,2);
        %hist(p_ymax(p_ymax<100))
        %title(['ymax ' num2str(ind_patient)])
        histogram(p_ymax(p_ymax<60),15)
        xlim([0,60])
    end
    
    %d
    figure(length(findobj('type','figure'))+3)
    for i=1:3
        subplot(3,1,i)
        p_d=params(m_pixels(:,4)==i,3);
        %title(['d ' num2str(ind_patient)])
        histogram(p_d(p_d<50),15)
        xlim([0,50])
    end
    
    %alpha
    figure(length(findobj('type','figure'))+4)
    for i=1:3
        subplot(3,1,i)
        p_alpha=params(m_pixels(:,4)==i,4);
        %title(['alpha ' num2str(ind_patient)])
        histogram(p_alpha(p_alpha<50),15)
        xlim([0,50])
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