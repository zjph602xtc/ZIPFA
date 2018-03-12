function [startvar] = EMzeropoisson_mat(data, tau, tol , maxiter, initialtau, Madj, display, initial, intercept)
% intercept: 'on' 'off'
% Madj: true false
% initial: [tau betas]

x = data(:,2:end);
y = data(:,1);
taubackup = tau;
if ~any(strcmp(initialtau,{'stable' 'initial' 'iteration'}))
    warning('Wrong initialtau, use initial mode instead.')
end

[n, nvar] = size(data);
if isempty(initial)
    if display
        disp('Initializing')
    end
    startvar = glmfit(x,y,'poisson','link','log','constant',intercept);
    startvar = [tau; startvar]';
else
    startvar = initial;
end

if strcmp(intercept,'on')
    Dx = [ones(n,1) x];
else
    Dx = x;
end

i = 1;
if display
    disp('Start maximizing ...')
end

while (i < maxiter && (i==1 || ...
        sqrt(sum((startvar(i,:)-startvar(i-1,:)).^2))/sqrt(sum(startvar(i,:).^2))>tol))
    % E step
    Z = zeros(1, n);
    Izero = y==0;
    if Madj
        m = evalin('caller','mi');
        Z(Izero) = 1./(1 + exp( startvar(i,1).*(Dx(Izero,:) * startvar(i,2:end)')-exp(m(Izero).*(Dx(Izero,:) * startvar(i,2:end)'))  ));
    else
        Z(Izero) = 1./(1 + exp( startvar(i,1).*(Dx(Izero,:) * startvar(i,2:end)')-exp(Dx(Izero,:) * startvar(i,2:end)')  ));
    end

    
    % M step
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,...
    'HessianFcn','objective','Display','off');
    if (i ~=1 && strcmp(initialtau,'stable'))
        tau = mean(log( n./ sum(Z==0) - 1) ./ (Dx * startvar(i,2:end)'));
    end
    if (strcmp(initialtau,'iteration'))
        tau = startvar(i, 1);
    end
  
    if Madj
        startvar = [startvar; fminunc(@lnL, [tau startvar(i,2:end)]', options, Dx,y,Z,m)'];
        %         disp('Madj')
    else
        startvar = [startvar; fminunc(@lnL, [tau startvar(i,2:end)]', options, Dx,y,Z,1)'];
%         startvar=[startvar; fsolve(@dlnL,[tau startvar(i,2:end)]', optimoptions('fsolve','display','off',...
%             'Algorithm', 'levenberg-marquardt','StepTolerance',1e-4,'InitDamping',50),Dx,y,Z)'];
    end
    
    i = i+1;
    if (strcmp('interation',initialtau) && abs(max(diff(startvar(:,1))))>50)
        disp('Try stable !')
        return %%continue
    end
    if display
        fprintf('This is %.0f th iteration, Frobenius norm diff = %g. \n',...
             i, sqrt(sum((startvar(i,:)-startvar(i-1,:)).^2))/sqrt(sum(startvar(i,:).^2)))
    end
end

end