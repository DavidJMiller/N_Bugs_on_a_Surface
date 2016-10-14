function [posFinal] = bugs_torus(options, prams)

%{
=================FUNCTIONS=================
function eqn: equation of the graph the particles are constrained to move on
function projTangent: projects the particles along the tangent 
function projSphere: projects the particles along the tangent back onto sphere
===========================================
%}

function [x,y,z] = eqn(theta, phi)
  x = (prams.Radius + prams.radius*cos(theta))*cos(phi);
  y = (prams.Radius + prams.radius*cos(theta))*sin(phi);
  z = prams.radius*sin(theta);
end

function [xCylinder,yCylinder,zCylinder] = cylinderEqn(theta, phi)
  xCylinder = prams.radius*cos(theta);
  yCylinder = prams.radius*sin(theta);
  zCylinder = mod(phi, 2*pi);
end 

function vel = projTangent(vector, theta, phi)
  x = -1*sin(phi);
  y = cos(phi);
  z = 0;
  X = -1*sin(theta).*cos(phi);
  Y = -1*sin(theta).*sin(phi);
  Z  = cos(theta);  
  t1 = [x; y; z];
  t2 = [X; Y; Z];
  normal = cross(t1, t2);
  vel = vector - normal*dot(vector, normal);
end

function [newTheta, newPhi] = projTorus(theta, phi, x, y, z)
  pair = [theta; phi];
  for l = 1:prams.newtonSteps
    theta = pair(1,1);
    phi = pair(2,1);
    dF_dTheta = 2*prams.radius*(sin(theta).*(-prams.Radius + x*cos(phi) + ...
      y*sin(phi)) - z*cos(theta));
    dF_dPhi = 2*(prams.radius*cos(theta) + prams.Radius).*(x*sin(phi) - y*cos(phi));
    d2F_dThetadPhi = 2*prams.radius*sin(theta).*(y*cos(phi) - x*sin(phi));
    d2F_d2Theta = 2*prams.radius*(cos(theta).*(-prams.Radius + x*cos(phi) + ...
      y*sin(phi)) + z*sin(theta));
    d2F_d2Phi = 2*(prams.radius*cos(theta) + prams.Radius).*(x*cos(phi) + y*sin(phi));
    F = [dF_dTheta; dF_dPhi];
    inverseJacobian = inv([d2F_d2Theta d2F_dThetadPhi; d2F_dThetadPhi d2F_d2Phi]);
    pair = pair - inverseJacobian*F;
  end
  newTheta = pair(1,1);
  newPhi = pair(2,1);
end

%Initialize the positions of the particles
%angles is a 2xN matrix  where the first row is the theta angle
%and the second row is the phi angle. Column m corresponds to particle m.
angles = zeros(2, prams.N);
for n = 1:prams.N 
  angles(1, n) = 2*pi*rand;
  angles(2, n) = 2*pi*rand;
  thetaInit = angles(1, n);
  phiInit = angles(2, n);
  [xInit,yInit,zInit] = eqn(thetaInit, phiInit);
  [xCylinderInit, yCylinderInit, zCylinderInit] = cylinderEqn(thetaInit, phiInit);
  posInit(1, n) = xInit;
  posInit(2, n) = yInit;
  posInit(3, n) = zInit; %+ floor(3*n/prams.N)*.1;
  posSquareInit(1, n) = thetaInit;
  posSquareInit(2, n) = phiInit;
  posCylinderInit(1, n) = xCylinderInit;
  posCylinderInit(2, n) = yCylinderInit;
  posCylinderInit(3, n) = zCylinderInit;
end

pos = posInit;
posSquare = posSquareInit;
posCylinder = posCylinderInit;

vel = zeros(3, prams.N);

for timeStep = 1:prams.ntime

  for k = 1:prams.N-1
  	direction = pos(:, k+1) - pos(:, k);
    vel(:, k) = projTangent(direction, angles(1,k), angles(2,k));
  end
  direction = pos(:, 1) - pos(:, prams.N);
  vel(:, prams.N) = projTangent(direction, angles(1,prams.N), angles(2,prams.N));

  if (options.record) && (timeStep == 1)
    writerObj = VideoWriter('Torus.avi'); % Name it.
    writerObj.FrameRate = 20; % How many frames per second.
    open(writerObj); 
  end

  %Plotting the simulation
  if options.iplot
    if timeStep == 1
      fId = figure;
    end 
    thetaGrid = linspace(0,2*pi);
    phiGrid = linspace(0,2*pi);
    [thetaGrid,phiGrid] = meshgrid(thetaGrid, phiGrid);
    subplot(2,11,[1 3])
    hold off;
    xAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*cos(phiGrid);
    yAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*sin(phiGrid);
    zAxis = prams.radius.*sin(thetaGrid);
    surf(xAxis,yAxis,zAxis, 'EdgeColor', 'none');
    hold on;
    plot3(pos(1,:), pos(2,:), pos(3,:), 'r.', 'markersize', 20);
    title('Bottom of Torus');
    hold on;
    axis equal
    axis([-(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.radius + .1) (prams.radius + .1)]);
    view(0, -90);
    subplot(2,11,[5 7])
    hold off;
    xAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*cos(phiGrid);
    yAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*sin(phiGrid);
    zAxis = prams.radius.*sin(thetaGrid);
    surf(xAxis,yAxis,zAxis, 'EdgeColor', 'none');
    hold on;
    plot3(pos(1,:), pos(2,:), pos(3,:), 'r.', 'markersize', 20);
    str = sprintf('N = %i' , prams.N);
    title(str);
    hold on;
    axis equal
    axis([-(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.radius + .1) (prams.radius + .1)]);
    subplot(2,11,[9 11])
    hold off;
    xAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*cos(phiGrid);
    yAxis = (prams.Radius + prams.radius*cos(thetaGrid)).*sin(phiGrid);
    zAxis = prams.radius.*sin(thetaGrid);
    surf(xAxis,yAxis,zAxis, 'EdgeColor', 'none');
    hold on;
    plot3(pos(1,:), pos(2,:), pos(3,:), 'r.', 'markersize', 20);
    title('Top of Torus');
    hold on;
    axis equal
    axis([-(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.Radius + prams.radius + .1) (prams.Radius + prams.radius + .1) ...
      -(prams.radius + .1) (prams.radius + .1)]);
    view(0, 90);
    subplot(2,11,[13 16])
    hold off;
    plot(posSquare(1,:), posSquare(2,:), 'r.', 'markersize', 20);
    hold on;
    axis equal 
    axis([0 2*pi 0 2*pi]);
    subplot(2,11,[18 21])
    hold off;
    xCylinderAxis = prams.radius*cos(thetaGrid);
    yCylinderAxis = prams.radius*sin(thetaGrid);
    zCylinderAxis = phiGrid;
    surf(xCylinderAxis, yCylinderAxis, zCylinderAxis, 'EdgeColor', 'none');
    hold on;
    plot3(posCylinder(1,:), posCylinder(2,:), posCylinder(3,:), 'r.', 'markersize', 20);
    hold on; 
    axis equal;
    axis([-prams.radius prams.radius -prams.radius prams.radius 0 2*pi]);
    pause(1e-2);
    if options.record
      frame = getframe(gcf);
      writeVideo(writerObj, frame);
    end
  end

  pos = pos + vel*prams.dt;
  
  for i = 1:prams.N
    thetaOld = angles(1,i);
    phiOld = angles(2,i);
    [angles(1,i), angles(2,i)] = ...
      projTorus(thetaOld, phiOld, pos(1,i), pos(2,i), pos(3,i));
  end
  
  for i = 1:prams.N
    [xNew,yNew,zNew] = eqn(angles(1,i), angles(2,i));
    pos(1,i) = xNew;
    pos(2,i) = yNew;
    pos(3,i) = zNew;
    [xCylinderNew, yCylinderNew,zCylinderNew] = cylinderEqn(angles(1,i), angles(2,i));
    posCylinder(1,i) = xCylinderNew;
    posCylinder(2,i) = yCylinderNew;
    posCylinder(3,i) = zCylinderNew;
    posSquare(1,i) = mod(angles(1,i), 2*pi);
    posSquare(2,i) = mod(angles(2,i), 2*pi);
  end
end

posFinal = pos;

if options.record
  close(writerObj);
end

clf;
close;

end 
