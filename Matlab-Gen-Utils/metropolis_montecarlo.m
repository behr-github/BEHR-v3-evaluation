function [ f, naccepts, naccepts_uphill, naccepts_downhill ] = metropolis_montecarlo( P, d, nsteps, f0, thin, burnin, lb, ub, Aineq, bineq, nonlincon )
%METROPOLIS_MONTECARLO Carries out metropolis monte carlo sampling on P
%   Detailed explanation goes here

nvar = numel(d);
f = nan(round((nsteps-burnin)/thin),nvar);

if ~exist('lb','var')
    lb = [];
end
if ~exist('ub','var')
    ub = [];
end

% Set up the initial starting point. If bounds are given, start anywhere in
% those bounds. If not, start somewhere assuming that the sample space
% should be centered on 0 and that the step size is ~100x smaller than the
% sample space. This is a weak assumption, which is why its much better if
% bounds can be given.
if ~exist('f0','var') || isempty(f0)
    if ~isempty(lb) && ~isempty(ub)
        fcurr = ((rand(1,nvar) - 0.5) .* (ub - lb)) + lb;
    else
        fcurr = (rand(1,nvar) - 0.5) .* 200 .* d;
    end
else
    fcurr = f0;
end

% By default assume no burnin and no thinning
if ~exist('burnin','var')
    burnin = 0;
end
if ~exist('thin','var')
    thin = 1;
end

% Make default bounds if necessary
if ~exist('lb','var') || isempty(lb)
    lb = -Inf(1,nvar);
end
if ~exist('ub','var') || isempty(ub)
    ub = Inf(1,nvar);
end
if ~exist('Aineq','var') || ~exist('bineq','var')
    Aineq = [];
    bineq = [];
end
if ~exist('nonlincon','var')
    nonlincon = [];
end

Pcurr = P(fcurr);
naccepts_downhill = 0;
naccepts_uphill = 0;
i=1;
for t=1:nsteps-1
    % First calculate the new possible coordinate.
    step = (rand(1,nvar) - 0.5) .* d .* 2;
    ftest = fcurr + step;
    
    % Decide whether to accept or reject the step
    Ptest = P(ftest);
    w = exp(-Ptest) / exp(-Pcurr);
    if any(ftest < lb) || any(ftest > ub)
        % Move not accepted: would exceed bounds
        rejectMove;
    elseif ~isempty(Aineq) && ~isempty(bineq) && any(Aineq * ftest) > bineq
        % Move not accepted: would excced linear constraints
        rejectMove;
    elseif ~isempty(nonlincon) && any(nonlincon > 0)
        % Move not accepted: would violate nonlinear constraints.
        rejectMove;
    %elseif Ptest < Pcurr;
    elseif exp(-(Ptest-Pcurr)) > rand(1)
        % "Downhill" move: accept automatically. Otherwise allow there to
        % be a random chance of an uphill move.
        naccepts_downhill = acceptMove(naccepts_downhill);
    elseif w <= rand(1)
        naccepts_uphill = acceptMove(naccepts_uphill);
    else
        % Move not accepted: stay here.
        rejectMove;
    end
end

naccepts = naccepts_downhill + naccepts_uphill;

    function nacc = acceptMove(nacc)
        if t > burnin && mod(t,thin) == 0
            f(i,:) = ftest;
        end
        nacc = nacc + 1;
        fcurr = ftest;
        Pcurr = Ptest;
        i=i+1;
    end

    function rejectMove()
        if t > burnin && mod(t,thin) == 0
            f(i,:) = fcurr;
        end
        i=i+1;
    end


end

