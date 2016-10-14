prams.N = 8;       % number of bugs
prams.T = 20;      % time horizon
prams.ntime = 200; % number of time steps

options.iplot = true; % plot the simulation with time

theta = (0:prams.N-1)'*2*pi/prams.N;
perturb = .5e0;
if (perturb > pi/prams.N)
  fprintf('WARNING: Bugs may not be in order: Pertubation is too large\n')
end

theta = theta + perturb*2*(rand(prams.N,1) - 1/2);
% randomly perturb the initial lcoation

zInit = exp(1i*theta);
% initial bug locations

[zFinal,min_dist]  = bugs_circle(options,prams,zInit);
% run bug simulation.  Return the final confiugration and the minimum
% distance between two succcessive bugs.  min_dist going to 0 indicates
% they've reached a fixed point.  min_dist going to
% sqrt(2-2*cos(2*pi/prams.N)) indicates that they are going towards a
% periodic orbit


periodic_dist = sqrt(2-2*cos(2*pi/prams.N));
fixed_dist = 0;
if 0.95*(periodic_dist - min_dist(end)) < (periodic_dist - min_dist(1))
  fprintf('Periodic orbit\n')
elseif (min_dist(end) - fixed_dist) < 0.95*(min_dist(1) - fixed_dist)
  fprintf('Fixed point\n')
else
  fprintf('Inconclusive\n')
end


