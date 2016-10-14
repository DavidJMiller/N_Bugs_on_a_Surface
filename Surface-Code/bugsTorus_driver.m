%{
=============PARAMETER SETTINGS=============
prams.NBegin: number of particles to start off with
prams.N: how many bugs you want to run the program up to
prams.perturbations: how many perturbations to iterate through
  range of perturbations is [0, π/prams.perturbations, 2π/prams.perturbations, ...,  π]
prams.lobeNumber: the number of lobes
prams.radius: the radius of the sphere
prams.runs: number of simulations for each permutation of prams.N & prams.perturbations
prams.dt: size of timestep
prams.ntime: number of time steps
prams.newtonSteps: number of iterations for Newton's Method
options.iplot: true if want to see simulation of particles moving, false othewise
============================================
%} 
prams.NBegin = 40;
prams.N = 2500;
prams.Radius = 2;
prams.radius = 1;
prams.perturbations = 10;
prams.runs = 1;
prams.dt = 0.1;
prams.ntime = 1500;
prams.newtonSteps = 3;
prams.delta = 1e-10;
options.iplot = true;
options.record = false;

if prams.Radius > prams.radius
  torusResults = zeros(prams.N-prams.NBegin+1, 1); 

  for N = prams.NBegin:prams.N
    prams.N = N;
    disp(N);
    for m = 1:prams.runs
      [posFinal] = bugs_torus(options, prams);
    end
      torusResults(prams.N-prams.NBegin+1,1) = total/prams.runs
  end
elseif prams.Radius == prams.radius
	fprintf('This is a Horn Torus! Program Terminated!\n');
else
	fprintf('This is a Spindle Torus! Program Terminated!\n');	
end
