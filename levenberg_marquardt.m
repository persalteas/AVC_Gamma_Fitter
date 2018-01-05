

function a=levenberg_marquardt(X, Y, a, lamb, nmax)
    nbfig=length(findobj('type','figure'));
	n = 1;
	r_lamb = lamb*ones(1,nmax);
	r_a = a*ones(4,nmax);
	r_grad = zeros(4,nmax);
	r_sse = zeros(1,nmax);
	while n<nmax
		grad = gcout(X, Y, a); % 1x4 vector
		r_grad(:,n) = grad;
		secderiv = hcout(X, a); % 4x4 matrix
		fprintf('============ n=%d ===========\na =',n)
        disp(a)
        nbfig = nbfig+1;
        figure(nbfig);
        hold on
        plot(X, Y, 'r')
        xlabel('Temps')
        ylabel('Concentration')
		while true
            % Calcul de la direction du mouvement
			HLM = (1+lamb)*secderiv;
			dk = - grad/HLM;
            ndk = norm(dk);
            fprintf('tentative lamb = %f  |dk| = %.10f  a+dk = ',lamb,ndk)
            disp(a+dk)
            
            % Calcul de la fonction gamma avant et après le mouvement
            la = zeros(size(X));
            ladk = zeros(size(X));
            for i=1:length(X)
                la(i) = gamma(X(i),a);
                ladk(i) = gamma(X(i), a+dk);
            end
            plot(X, la, 'b', X, ladk, 'g');
            
			if ndk>20
				disp('explosion !')
				break
			end
			if ndk<1e-9
				disp('eteint...')
				break
			end
			if cout(X, Y, a+dk) < cout(X, Y, a)
				a = a+dk;
				r_a(:,n) = a;
				lamb = lamb/10;
				r_lamb(n) = lamb;
				break
			else
				lamb = lamb*10;
				r_lamb(n)=lamb;
			end
		end
		r_sse(n) = cout(X, Y, a);
		if (ndk>20 ||  ndk<1e-9)
				break
        end
        hold off;
		n = n + 1;
	end
	if n==nmax
		disp('poh trouve...');
	end
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

function f=cout(X,Y,a)
    distances = zeros(size(X));
    for i=1:size(X)
        distances = (Y(i) - gamma(X(i),a))^2;
    end
    f = 0.5*sum(distances);
end

function f=gcout(X,Y, a)
    tmax = a(1);
    ymax = a(2);
    d = a(3);
    alpha = a(4);
    f = [0 0 0 0];
    for i=1:length(X)
        t = X(i);
        y = Y(i);
        if (t<=d)
            v = [0 0 0 0];
        else
            T = (t-d)/tmax;
            dfdd = -ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax;
            dfdalpha = ymax*exp(alpha*(1-T))*T^alpha*(log(T)+1-T);
            dfdymax = T^alpha*exp(alpha*(1-T));
            dfdtmax = d*ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax/tmax;
            distance = (y - gamma(t,a));
            v = [distance*dfdtmax distance*dfdymax distance*dfdd distance*dfdalpha];
        end
        f = f+v;
    end
    f = -f;
end

function f=hcout(X,a)
    tmax = a(1);
    ymax = a(2);
    d = a(3);
    alpha = a(4);
    f = zeros(4,4);
    for i=1:length(X)
        t = X(i);
        if (t<=d)
            v = zeros(4,4);
        else
            T = (t-d)/tmax;
            dfdtmax = d*ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax/tmax;
            dfdymax = T^alpha*exp(alpha*(1-T));
            dfdd = -ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax;
            dfdalpha = ymax*exp(alpha*(1-T))*T^alpha*(log(T)+1-T);
            v = [ dfdtmax*dfdtmax dfdtmax*dfdymax dfdtmax*dfdd dfdtmax*dfdalpha;
                  dfdymax*dfdtmax dfdymax*dfdymax dfdymax*dfdd dfdymax*dfdalpha;
                  dfdd*dfdtmax dfdd*dfdymax dfdd*dfdd dfdd*dfdalpha;
                  dfdalpha*dfdtmax dfdalpha*dfdymax dfdalpha*dfdd dfdalpha*dfdalpha ];
        end
        f = f+v;
    end
end


% function f=ggamma(t, a)
%     tmax = a(1);
%     ymax = a(2);
%     d = a(3);
%     alpha = a(4);
%     if (t<=d)
%         f = 0;
%     else
%         T = (t-d)/tmax;
%         f = ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax;
%     end
% end

% function f=hgamma(ymax, tmax, alpha, d, t)
%     if (t<=d)
%         f = 0;
%     else
%         T = (t-d)/tmax;
%         f = ymax*alpha*T^(alpha-2)*exp(alpha*(1-T))*(alpha-1-2*T*alpha+alpha*T*T)/tmax/tmax;
%     end
% end

