

function a=levenberg_marquardt(X, Y, a, lamb, nmax)
	n = 1;
	while n<nmax
        grad = gcout(X, Y, a); % 4x1 vector
        secderiv = hcout(X, a); % 4x4 matrix
        disp(['====== n=',num2str(n)]);
        while (1)
			HLM = (ones(4)+diag(lamb*ones(1,4))).*secderiv;
			dk = - grad/HLM;
			disp(['l = ',num2str(lamb)])
            disp(['a = ', num2str(a)])
            disp(['DLM = ', num2str(dk)])
%             if norm(dk)>20
% 				disp('explosion !')
% 				break
%             end
            if norm(dk)<0.0000001
                disp('eteint...')
                break
            end
            if a(1)+dk(1)>0
                cout(X, Y, a)
                cout(X, Y, a+dk )
                if cout(X, Y, a+dk ) < cout(X, Y, a)
                    a = a+dk;
                    lamb = lamb/10;
                    break
                else
                    lamb = lamb*10;
                end
            else
                lamb=lamb*10;
            end
            
        end
        %if (norm(dk)>20 || norm(dk)<0.000000001)
        if (norm(dk)<0.0000001)
				break
        end
        n = n + 1;
    end
	if n==nmax
		disp('poh trouve...');
	end
    
end

function f=lambda(t, a)
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
    for i=1:length(X)
        distances(i) = (Y(i) - lambda(X(i),a))^2;
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
            dfdtmax = ymax*alpha*T^(alpha)*exp(alpha*(1-T))*(T-1)/tmax;
            dfdymax = T^alpha*exp(alpha*(1-T));
            dfdd = -ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax;
            dfdalpha = ymax*exp(alpha*(1-T))*T^alpha*(log(T)+1-T);
            distance = (y - lambda(t,a));
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
            %dfdtmax = d*ymax*alpha*T^(alpha-1)*exp(alpha*(1-T))*(1-T)/tmax/tmax;
            dfdtmax = ymax*alpha*T^(alpha)*exp(alpha*(1-T))*(T-1)/tmax;
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

% function f=glambda(t, a)
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

% function f=hlambda(ymax, tmax, alpha, d, t)
%     if (t<=d)
%         f = 0;
%     else
%         T = (t-d)/tmax;
%         f = ymax*alpha*T^(alpha-2)*exp(alpha*(1-T))*(alpha-1-2*T*alpha+alpha*T*T)/tmax/tmax;
%     end
% end

