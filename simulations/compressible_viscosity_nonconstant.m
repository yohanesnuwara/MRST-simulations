############################################################
#    Compressible Simulation on 3D Homogeneous Reservoir   #
#            Pressure-dependent Viscosity - Gas            #
#                                                          #
# Compare to the constant visco. Lines that are different: #
# * Visco definition in line 26 and 27                     #
# * Darcy equation in line 97                              #
# * Flow rate through each well connection in line 115     # 
############################################################

## Setup model and solve the initial hydrostatic pressure

# Define geometry
[nx,ny,nz] = deal( 10, 10, 10);
[Lx,Ly,Lz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

# Define rock
rock = makeRock(G, 30*milli*darcy, 0.3);
%mrstModule add spe10
%rock = getSPE10rock(41:50,101:110,1:10);

# Rock compressibility is constant, pore volume is function of pressure (Eq 7.6)
cr = 1e-6/barsa;
p_r = 200*barsa;
pv_r = poreVolume(G, rock);
pv = @(p) pv_r .* exp( cr * (p - p_r) );

# Compressibility is constant, viscosity is pressure dependent 
# rho function of pressure (Eq 7.7)
[mu0,c_mu] = deal(5*centi*poise, 2e-3/barsa);
mu = @(p) mu0*(1+c_mu*(p-p_r));

c = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS = 750*kilogram/meter^3;
rho = @(p) rho_r .* exp( c * (p - p_r) );

# Define horizontal well and its perforations
nperf = 8;
pwf = 100*barsa; % BHP of well
I = repmat(2, [nperf, 1]);
J = (1:nperf).'+1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'producer', 'Dir', 'x');

# solve initial hydrostatic pressure
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);

%% NOTE: deval function is not available in Octave. However, I have added
%% my own function, originally written by Andres Codas in MRST forum
mrstModule add nuwara

p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);

# Plot initial hydrostatic pressure
figure(1)
show = true(G.cells.num,1);
cellInx = sub2ind(G.cartDims, ...
[I-1; I-1; I; I; I(1:2)-1], ...
[J ; J; J; J; nperf+[2;2]], ...
[K-1; K; K; K-1; K(1:2)-[0; 1]]);
show(cellInx) = false;
plotCellData(G,p_init/barsa, show, 'EdgeColor','k'), colorbar;
title('Initial Hydrostatic Pressure [bar]');
colormap (jet(55))
plotWell(G,W, 'height',10);
view(-125,20);

## Discretize model

# Compute map between interior faces and cells
C = double(G.faces.neighbors);
%% NOTE: In the book p. 207, intInx is not defined. Here, I use another source
%% in website and find the 2014 MRST book documented this. So, I modified the code.
intInx = all(C ~= 0, 2);
C = C(intInx, :);

%% NOTE: exterior faces are not included because assumed no flow

# Define averaging operator
n = size(C,1);
D = sparse([(1:n)'; (1:n)'], C, ...
ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x) D*x;
div = @(x) -D'*x;
avg = @(x) 0.5 * (x(C(:,1)) + x(C(:,2)));

# Compute transmissibilities
hT = computeTrans(G, rock); % Half-transmissibilities
cf = G.cells.faces(:,1);
nf = G.faces.num;
T = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T = T(intInx); % Restricted to interior

# Darcy's equation
gradz = grad(G.cells.centroids(:,3));
v = @(p) -(T./mu(avg(p))).*( grad(p) - g*avg(rho(p)).*gradz );

# Continuity equation for each cell C
presEq = @(p,p0,dt) (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) ...
+ div( avg(rho(p)).*v(p) );

%% NOTE: until this point, when I compute the norm, result is different from book
%% p. 208. My result 8.3e-4. Book result is 1.5e-6. I hope it's right.
display(norm(v(p_init))*day);

## Define well model
%% NOTE: pay attention to statement in p. 208 about well model
%% The following code is hence general for all well geometries (vert, hor, or dev)

wc = W(1).cells; % connection grid cells
WI = W(1).WI; % well-indices
dz = W(1).dZ; % depth relative to bottom-hole
p_conn = @(bhp) bhp + g*dz.*rho(bhp); %connection pressures
q_conn = q_conn = @(p,bhp) WI.*(rho(p(wc))./ mu(p(wc))) .* (p_conn(bhp) - p(wc));

# Compute total volumetric well rate
rateEq = @(p, bhp, qS) qS-sum(q_conn(p, bhp))/rhoS;

# Declare the well condition as constant BHP
ctrlEq = @(bhp) bhp-pwf;

## Initialize simulation loop

# Initialize AD variables
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);

# Set indices
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
nc = G.cells.num;
[pIx, bhpIx, qSIx] = deal(1:nc, nc+1, nc+2);

# Set timesteps
[numSteps, totTime] = deal(52, 365*day); % time-steps/ total simulation time
[tol, maxits] = deal(1e-5, 10); % Newton tolerance / maximum Newton its
dt = totTime / numSteps;

sol = repmat(struct('time',[], 'pressure',[], 'bhp',[], 'qS',[]), [numSteps+1, 1]);
sol(1) = struct('time', 0, 'pressure', value(p_ad), ...
'bhp', value(bhp_ad), 'qS', value(qS_ad));

## Main simulation loop

t = 0; step = 0;
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));
   % Newton loop
   resNorm = 1e99;
   p0  = value(p_ad); % Previous step pressure
   nit = 0;
   while (resNorm > tol) && (nit <= maxits)
      % Add source terms to homogeneous pressure equation:
      eq1     = presEq(p_ad, p0, dt);
      eq1(wc) = eq1(wc) - q_conn(p_ad, bhp_ad);
      % Collect all equations
      eqs = {eq1, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      bhp_ad.val = bhp_ad.val + upd(bhpIx);
      qS_ad.val  = qS_ad.val  + upd(qSIx);

      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
                            'bhp', value(bhp_ad), 'qS', value(qS_ad));
   end
end

# Plot production rate and pressure decrease
addpath 'C:\Octave\Octave-5.2.0\mingw64\share\octave\5.2.0\m\plot\draw'

figure(2)
[ha,hr,hp] = plotyy(...
   [sol(2:end).time]/day, -[sol(2:end).qS]*day, ...
   [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), 'stairs', 'plot');
set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
set(ha(2),'YLim',[100 210],'YTick',100:50:200);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
ylabel(ha(2), 'avg pressure [bar]');

# Plot evolution of pressure distribution
%% NOTE: Ignore the warning. It's related to the well plot, but it's alright.
%clf;
figure(3)
steps = [2 5 10 20];
for i=1:4
   subplot(2,2,i);
   set(gca,'Clipping','off');
   plotCellData(G, sol(steps(i)).pressure/barsa, show,'EdgeColor',.5*[1 1 1]);
   plotWell(G,W);
   view(-125,20)
   %caxis([115 205]);
   axis tight off, colorbar('southoutside');
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',14);
end

colormap(jet(55));