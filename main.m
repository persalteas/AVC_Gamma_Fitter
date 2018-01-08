%% main 
clc;
warning('off', 'MATLAB:nearlySingularMatrix'); % Ne pas afficher les warnings de matrice presque singulière: osef
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart'); % Ne pas afficher les warnings sur les nombres complexes.
%addpath('NIfTI_20140122');


% display = input('Afficher des sorties ? \n0 = non | 1 = texte | 2 = graphes | 3 = tout \n:');
% if (display >=2), disp('Pour passer au pixel suivant, cliquez sur le graphe.'), end
% patient = input('\nNumero du patient ? (3 caracteres) : ','s');
% AnalysePatient(patient, display);

%AnalysePatient('048',1);
AnalysePatient('249',1);
%AnalysePatient('170',1);
%AnalysePatient('260',1);

%% A DISCUTER LA PROCHAINE FOIS:
%
%	- Pourquoi d = find(mean_1>=0.4,1) comme valeur initiale ?
%	- Pourquoi la fonction de cout c'est 1/2 de la somme des écarts (et pas	1/N) ?
%	- Pourquoi while any(secderiv(:)) ?
%	- Le projet demande d'utiliser R² comme fonction de cout non ?
%		Nous on l'a fait en minimisant 0.5*sum(Yt - f(t))
%		J'ai l'impression qu'elle veut qu'on minimise sum(Yt-<Yt>)/sum(Yt-f(t))
%		qui est pourtant indéfinie quand Yt=f(t), notamment quand les deux sont nulles
%
