function [ repby_j, normentro, turntimes, rectimes, brrectime,...
    entropycols, dobcoeff, damp, convergerate, halflife, nlogdobcoeff,...
    w1] = calc_communityProps( P, e, P_interval )
%CALC_COMMUNITYPROPS Calculate successional indices (Hill et al. 2004)
%
%INPUT
%   P_samples - transition probabilities; row - sample, col - param
%   e - state number for bare space
%   P_interval - scale of P in months (ex. 12 if annual)
%
%OUTPUT
%   repby_j - average p of j replacing all other other states
%   normentro - normalized entropy
%   turntimes - vector of turnover times (inverse turnover rates); in years
%   rectimes - vector of recurrence times; in years
%   brrectime - recurrence time of empty state (CCA + turf)
%   entropycols - vector, entropy of functional group transitions
%   dobcoeff - Dobrushin's ergodicity coefficient
%   damp - damping ratio
%   convergerate - convergence rate 
%   halflife - half life  of x months
%   nlogdobcoeff - short-term convergence rate 
%   w1 - dominant right eigenvector

[l,W,V] = sorted_eigs(P);
s = size(P,1);
sbio = [1:s];
sbio(e) = [];  %take out the index for bare space
w1 = W(:,1);

% % % Disturbance
% % % disturb - p that an occupied point, randomly selected from the 
% % % stationary community, is disturbed between t and t+1
% disturb = 0;
% for i=1:(s-1)
%     state = sbio(i);
%     disturb = disturb + w1(state)*P(e,state);
% end
% disturb = disturb / sum(w1(sbio));

% % % Colonization
% % % p that an empty point is colonized by some func. group
% colon = 1-P(e,e);

% % % Persistence
% % % mean persistence rate
% meanpersist = 0;
% for i=1:s
%    meanpersist = meanpersist + w1(i)*P(i,i);
% end

% % % Replacement by, averaged over all states
% % % p of replacement by functional group i, averaged over all states
% repby = 0;
% for i = 1:s-1
%     state = sbio(i);
%     pijs = 0;
%     for j = 1:s-1
%         if ne(i,sbio(j))
%             pijs = pijs + P(i,sbio(j));
%         end
%     end
%     repby = repby + w1(state)*pijs;
% end
% repby = (1/(s-2)) * repby / sum(w1(sbio));

% Replacement by group j
repby_j = zeros(s,1);
for i = 1:s-1
    state = sbio(i);    
    repby_j(state) = (sum(P(state,sbio)) - P(state,state)) / (s - 2);
end

% % % Replacement of group j
% repof_j = zeros(s,1);
% for i = 1:s-1
%     state = sbio(i);
%     repof_j(state) = sum(P(sbio,state)) - P(state,state);
% end

% % % Replacement of, averaged over all biotic states
% % % p that a point, randomly selected from the stationary dist, 
% % % is replaced by a different func. group
% repof = 0;
% for i=1:s-1
%     state = sbio(i);
%     repof = repof + w1(state)*(1 - P(state,state) - P(e,state));
% end
% repof = repof / sum(w1(sbio));

% % % Proportion of bare space in stationary distribution
% bare = w1(e);

% Normalized entropy
normentro = 0;
for i=1:s
    entropy = 0;
    for j=1:s
        tmp = P(j,i)*log(P(j,i));
        if isnan(tmp)
            tmp = 0;
        end
        entropy = entropy + tmp;
    end
    normentro = normentro + w1(i)*entropy;
end
normentro = - normentro / -log(1/s);

% Turnover rates (vector)
turnrates = 1 - diag(P);

% Turnover times (vector)
turntimes = (1 ./ turnrates) * (P_interval / 12);   % rescaled to years

% % % Biotic turnover rate - averaged over all biotic states
% bioturnrate = w1(sbio) .* turnrates(sbio);
% bioturnrate = sum(bioturnrate) / sum(w1(sbio));

% % % Biotic turnover time - averaged over all biotic states
% bioturntime =  w1(sbio) .* turntimes(sbio);
% bioturntime = (sum(bioturntime) / sum(w1(sbio)));

% Recurrence times (vector)
rectimes = (1 - w1) ./ (w1 .* (1 - diag(P)));   % given in the units of matrix resolution
rectimes = rectimes * (P_interval / 12);        % rescaled

% % % Mean biotic recurrence time
% % % weight recurrence times according to stationary abundances
% mubiorectime = w1(sbio) .* rectimes(sbio);     
% mubiorectime = sum(mubiorectime) / sum(w1(sbio));

% % % Bare rock turnover rate
% brturnrate = turnrates(e);

% % %Bare rock turnover time (years)
% brturntime = turntimes(e);

%Bare rock recurrence time / mean time between disturbances 
brrectime = rectimes(e);

% % % Functional group turnovers between disturbances
% spturns = brrectime / bioturntime;

% Entropy of functional group transitions given group (vector)
% Inverse measure of predictability (High entropy, low predictability)
entropycols = zeros(s,1);
for i = 1:s
    tmp = P(:,i) .* log(P(:,i));
    tmp(find(P(:,i) == 0)) = 0;
    entropycols(i) = - sum(tmp) / -log(1/s);
end

% CONVERGENCE MEASURES
% Damping ratio
damp = l(1) / abs(l(2));  %rho

% Dobrushin's ergodicity coefficient
coeffs = [];
for i = 1:(s-1)
    for j = (i+1):s
        pdiff = P(:,i) - P(:,j);
        coeffs(end+1) = norm(pdiff,1);        
    end
end
dobcoeff = max(coeffs) / 2;

% convergence rate    x% per month  %22% per month? 
convergerate = log(damp);  

% half life  of x months
halflife = log(2) / log(damp);

% short-term convergence rate 
nlogdobcoeff = - log(dobcoeff);

% % % Biotic evenness - Shannon-Wiener diversity  
% Hwp = - sum(w1(sbio) .* log(w1(sbio)));
% J = Hwp / log(s-1);

% % %Proportional change in evenness
% dJ = zeros(s,1);
% for i=1:s
%     Pi = P; Ji = 0;
%     Pi(:,i) = [];
%     Pi(i,:) = [];
%     normalizer = repmat(sum(Pi),(s-1),1);
%     Pi = Pi ./ normalizer;
%     [li,wi,vi] = sorted_eigs(Pi);
%     wi = wi(:,1);
%     if ne(i,e)
%         if gt(i,e)
%             adj_e = e;
%         else
%             adj_e = e-1;
%         end
%         wi(adj_e) = [];     %take out empty state
%         Ji = - sum(wi .* log(wi)) / log(s-2);
%     else
%         Ji = - sum(wi .* log(wi)) / log(s-1);
%     end
%     dJ(i) = (Ji - J) / J;
% end
%======================================== 
    
    
% % disp(['Disturbance:  ', num2str(disturb)])
% % disp(['Colonization:  ', num2str( colon )])
% % disp(['Persistence:  ', num2str( meanpersist )]) 
% % disp(['Replacement by:  ',num2str( repby )]) 
% % disp(['Replacement of:  ',num2str( repof )]) 
% % disp(['Bare space:  ',num2str( bare )]) 
% % disp(['Biotic turnover rate:  ',num2str( bioturnrate )]) 
% % disp(['Biotic turnover time (years):  ',num2str( bioturntime )]) 
% % disp(['Mean biotic recurrence time (years):  ',num2str( mubiorectime )])
% % disp(['BR turnover rate:  ',num2str( brturnrate )]) 
% % disp(['BR turnover time (years):  ',num2str( brturntime )])  
% % disp(['Species turnovers:  ',num2str( spturns )])
% % disp(['Biotic evenness:  ',num2str(J)])
% % disp(['Proportional change in biotic evennness'])
% % disp([dJ]) %% uncomment

% disp(['Replacement by group j:  '])
% disp([repby_j])
% disp(['Normalized entropy:  ',num2str( normentro )]) 
% disp(['Turnover times (y)'])
% disp(turntimes)   %% uncomment
% disp(['Recurrence time (years)'])
% disp([rectimes])   %% uncomment
% disp(['BR recurrence time (years):  ',num2str( brrectime )])
% disp(['Entropy of transitions'])
% disp([entropycols])   %% uncomment
% 
% disp(['Damping ratio:  ',num2str(damp)]) 
% disp(['convergence rate: log(rho) = ',num2str(convergerate),' x100 percent per month'])
% disp(['half-life: log(2)/log(rho) = ',num2str(halflife),' months'])
% disp(['Dobrushins coefficient (alpha):  ',num2str( dobcoeff )]) 
% disp(['Short-term convergence rate: -log(alpha) = ',num2str(nlogdobcoeff),' x100 percent per month'])
% 



end

