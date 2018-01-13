%% main 
clc;
warning('off', 'MATLAB:nearlySingularMatrix'); % Ne pas afficher les warnings de matrice presque singulière: osef
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart'); % Ne pas afficher les warnings sur les nombres complexes.
%addpath('NIfTI_20140122');


% display = input('Afficher des sorties ? \n0 = non | 1 = texte | 2 = graphes | 3 = tout \n:');
% if (display >=2), disp('Pour passer au pixel suivant, cliquez sur le graphe.'), end
% patient = input('\nNumero du patient ? (3 caracteres) : ','s');
% AnalysePatient(patient, display);

AnalysePatient('048',1);
%AnalysePatient('249',1);
%AnalysePatient('170',1);
%AnalysePatient('260',3);
