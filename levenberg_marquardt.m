
%% Optimisation d'une fonction de cout par la methode de Levenerg-Marquardt:
% -------------------------------------------------------------------------
% FONCTION DE COUT A MINIMISER: la difference entre Y (observe) et gamma(X, params) (predit).

% La fonction de cout depend de 4 parametres: tmax, ymax, d et alpha.
% Elle prend donc ses antecedents dans un espace de 4 dimensions.
% La methode choisit des mouvements dans cet espace 4D e partir d'un point de depart
% de faeon e descendre les pentes de la fonction de cout, jusqu'e arriver dans un minimum.
%
% La pente dans les 4 directions de l'espace est donnee par le gradient gcout (derivees 1eres)
% ou la hessienne hcout (derivees 2ndes). En pratique, on construit une matrice HLM intermediaire entre
% les deux, en fonction d'un parametre lambda qui donne plus ou moins de poids au degre 1 ou 2.

function params=levenberg_marquardt(X, Y, params, lamb, nmax, numpixel, display)
	% X les valeurs d'abscisses
	% Y les valeurs d'intensite observees (ordonnees)
	% params les valeurs initiales des parametres
	% lambda valeur initiale de lambda (qui orient la methode L-M vers le degre 1 ou 2)
	% nmax nombre maximum d'iterations
	
	% Préparation d'un graphique
	if display>=2
		nbfig = length(findobj('type','figure')); % recupere le nombre de figures deje ouvertes par MatLab
		nbfig = nbfig+1;
		figure(nbfig);
		hold on
		ylim([ 0 max(Y)+10 ]);
		xlabel('Temps')
		ylabel('Concentration')
		title(sprintf('Optimisation de la distribution du pixel %d', numpixel));
	end
	
	subiteration_counter = 0;
	motif_arret = '';
	stop = 0;
	n = 1;
	while n<nmax % tant qu'on a le droit de faire des iterations
		
		% On calcule les derivees partielles 1eres et 2ndes 
		% de la fonction de cout avec les valeurs de parametres actuelles
		grad = gcout(X, Y, params); % 4x1 vector
		secderiv = hcout(X, params); % 4x4 matrix

		% while any(secderiv(:)) % pourquoi ?
		while true 
			subiteration_counter = subiteration_counter + 1;
			
			% calcul de la direction du mouvement dk:
			HLM = (ones(4)+diag(lamb*ones(1,4))) .* secderiv;
			dk = - grad/HLM;
			ndk = norm(dk);
			
			% les nouveaux parametres a essayer
            newparams = params + dk;

			% Calcul et affichage de la fonction gamma apres le mouvement propose
			if display >=2
				gamma_apres = zeros(size(X));
				for i=1:length(X)
					gamma_apres(i) = gamma(X(i),newparams);
				end
				plot(X, gamma_apres, 'g');
			end
				
			% On n'autorise que les mouvements a parametres positifs
			if prod(newparams>0)
				
				% Tests de la convergence
				if ndk < 10^-7 % Le mouvement propose est trop petit: l'algorithme s'est eteint, on s'arrete ici.
					motif_arret = 'extinction';
					stop = 1;
					break
				end
				if ndk > 200 % Le mouvement propose est gigantesque: il y a un probleme, on s'arrete
					motif_arret = 'explosion';
					stop = 1;
					break
				end
				if abs( cout(X,Y,newparams) - cout(X,Y,params) ) < 10^-7 % La difference de cout est faible : on s'arrete
					motif_arret = 'stagnation';
					stop = 1;
					break
				end
				
				% Accepter ou non les nouveaux parametres ?
				if cout(X,Y,newparams) < cout(X,Y,params) % le mouvement propose a bien permis de diminuer le cout
					params = params+dk;
					lamb = lamb/10; % ca fonctionne bien, on peut diminuer lambda (donner plus de poids au 1er degre)
					% disp('ca marche, go iteration suivante')
					break
				else
					lamb = lamb*10; % ca n'a pas marche, il faut faire des mouvements plus fins (donner plus de poids au 2nd degre)
				end
				
			else
				lamb = lamb*10; % comme si on refusait les parametres proposes.
			end
		end % iteration suivante.
		
		if (stop) % Sortir aussi du 2e while si on s'arrete
			break
		end
		n = n + 1;
	end
	if n == nmax % On a itere tres longtemps sans trouver de minimum
		if (display==1 || display==3), fprintf('pas de loi optimum trouvee en %d iterations !\t', n), end
	else % On s'arrete sur un truc pas mal
		if (display==1 || display==3), fprintf('pixel %d optimise en %d iterations (%d sous-iterations) par %s,\t', numpixel, n, subiteration_counter, motif_arret), end
	end
	if (display==1 || display==3), fprintf(' tmax = %6.2f\tymax = %6.2f\td = %6.2f\talpha = %6.2f \n', params(1), params(2), params(3), params(4)), end
	if display >= 2
		plot(X, Y, 'r')
		plot(X, gamma_apres, 'b');
		hold off;
		waitforbuttonpress; % Attend le clic de l'utilisateur sur le graphique pour passer au suivant
		close;
	end
end

% La fonction de la ditribution gamma
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

% La fonction de cout a optimiser: la moyenne des
% differences entre Y et gamma(X) au carre
function f=cout(X,Y,params)
	distances = zeros(size(X));
	for i=1:length(X)
		distances(i) = (Y(i) - gamma(X(i),params))^2;
	end
	f = 0.5*sum(distances); % pourquoi 0.5 ?
end

% Le gradient de la fonction de cout,
% Soit les derivees partielles (premieres) par rapport e 
% tmax, ymax, d et alpha
function f=gcout(X,Y,params)
	tmax = params(1);
	ymax = params(2);
	d = params(3);
	alpha = params(4);

	f = [0 0 0 0];
	for i=1:length(X)
		t = X(i);
		y = Y(i);
		if (t<=d) % euh ?
			v = [0 0 0 0];
		else
			T = (t-d)/tmax;
			dfdtmax = ymax * alpha * T^(alpha) * exp(alpha*(1-T)) * (T-1) / tmax;
			dfdymax = T^alpha * exp(alpha*(1-T));
			dfdd = - ymax * alpha * T^(alpha-1) * exp(alpha*(1-T)) * (1-T) / tmax;
			dfdalpha = ymax * exp(alpha*(1-T)) * T^alpha * (log(T)+1-T);
			distance = (y - gamma(t,params));
			v = [distance*dfdtmax distance*dfdymax distance*dfdd distance*dfdalpha];
		end
		f = f+v;
	end
	f = -f;
end

% La hessienne de la fonction de cout,
% Soit les derivees partielles secondes par rapport a
% tmax, ymax, d et alpha
function f=hcout(X,params)
	tmax = params(1);
	ymax = params(2);
	d = params(3);
	alpha = params(4);
	f = zeros(4,4);
	for i=1:length(X)
		t = X(i);
		if (t<=d)
			v = zeros(4,4);
		else
			T = (t-d)/tmax;
			dfdtmax = ymax * alpha * T^(alpha) * exp(alpha*(1-T)) * (T-1) / tmax;
			dfdymax = T^alpha * exp(alpha*(1-T));
			dfdd = - ymax * alpha * T^(alpha-1) * exp(alpha*(1-T)) * (1-T) / tmax;
			dfdalpha = ymax * exp(alpha*(1-T)) * T^alpha * (log(T)+1-T);
			v = [ dfdtmax*dfdtmax   dfdtmax*dfdymax	 dfdtmax*dfdd	dfdtmax*dfdalpha;
				  dfdymax*dfdtmax   dfdymax*dfdymax	 dfdymax*dfdd	dfdymax*dfdalpha;
				  dfdd*dfdtmax	  dfdd*dfdymax		dfdd*dfdd	   dfdd*dfdalpha;
				  dfdalpha*dfdtmax  dfdalpha*dfdymax	dfdalpha*dfdd   dfdalpha*dfdalpha ];
		end
		f = f+v;
    end
end

