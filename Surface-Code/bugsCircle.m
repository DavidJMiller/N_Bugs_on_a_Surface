function [zFinal,min_dist] = bugs_circle(options,prams,zInit)
% track bugs on a circle that are chasing each other.
% Everything can be (and is) done with complex variables since we are
% working with the unit circle.

dt = prams.T/prams.ntime;
% time step size

z = zInit;
% assign initial condition to z

min_dist = zeros(prams.ntime+1,1);
min_dist(1) = min(min(abs(z(2:end) - z(1:end-1))),...
    abs(z(1)-z(prams.N)));

vel = zeros(prams.N,1);

for j = 2:prams.ntime+1
  for k = 1:prams.N-1
    dir = z(k+1) - z(k);
    % direction between successive bugs (except betwen the last and the
    % first)

    if abs(dir) ~= 2
      vel(k) = -1i*imag(z(k)*conj(dir))*z(k);
    else
      vel(k) = 0;
    end
    % project the velocity onto the tangent plane of each bug

  end

  dir = z(1) - z(prams.N);
  % compute direction of last bug

  if abs(dir) ~= 2
    vel(prams.N) = -1i*imag(z(prams.N)*conj(dir))*z(prams.N);;
  else
    vel(prams.N) = 0;
  end

  vel(prams.N) = -1i*imag(z(prams.N)*conj(dir))*z(prams.N);
  % project the velocity of the last bug onto its tangent plane

  if options.iplot
    hold off;
    plot(exp(1i*(0:2*pi/50:2*pi)),'k');
    hold on;
    plot(z(:,1),'r.','markersize',20);
    axis equal
    pause(1e-2);
  end

  z = z + dt*vel;
  % move the bug along the tangent plane
  z = z./abs(z);
  % project the bug back onto the unit circle through the normal vector

  min_dist(j) = min(min(abs(z(2:end) - z(1:end-1))),...
      abs(z(1)-z(prams.N)));
  % compute minimum distance between successive bugs
end
zFinal = z;
