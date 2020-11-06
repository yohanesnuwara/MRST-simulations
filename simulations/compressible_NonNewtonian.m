############################################################
#    Compressible Simulation on 3D Homogeneous Reservoir   #
#              Non-Newtonial Fluid - Polymer               #
#                                                          #
#    Two methods: Cell-based, or Face-based. "switch_to="  #
############################################################

%% NOTE: Face-based is more robust and speedy than cell-based as it solves the 
%% problems in cell-based method outlined in p. 217 last paragraph. However, 
%% its potential drawback is numerical smearing. 

mrstModule add nuwara

# Choose method
switch_to='face_based'; % 'cell_based' or 'face_based'

# Inputs
[nx,ny,nz] = deal( 10,  10, 10);
[Lx,Ly,Lz] = deal(200, 200, 50);

rock = makeRock(G, 30*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;

c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;

numSteps = 52;
totTime  = 365*day;
tol      = 1e-5;
maxits   = 100;

%% Well model
nperf = 8;
pwf = 300*barsa; % BHP of well
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', 'x' );


if switch_to=='cell_based'
  %% Define geometric quantitites
  % Grid that represents the reservoir geometry
  G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
  G = computeGeometry(G);

  % Discrete operators
  N  = double(G.faces.neighbors);
  intInx = all(N ~= 0, 2);
  N  = N(intInx, :);
  n = size(N,1);
  C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
  aC = bsxfun(@rdivide,0.5*abs(C),G.faces.areas(intInx))';

  grad = @(x) C*x;
  div  = @(x) -C'*x;
  cavg = @(x) aC*x;
  favg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
  clear aC C N;

  %% Rock model and transmissibilities
  pv_r = poreVolume(G, rock);
  pv   = @(p) pv_r .* exp( cr * (p - p_r) );
  clear pv_r;

  hT = computeTrans(G, rock);
  cf = G.cells.faces(:,1);
  nf = G.faces.num;
  T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
  T  = T(intInx);
  clear hT;

  %% Basic fluid model
  rho   = @(p) rho_r .* exp( c * (p - p_r) );
  if exist('fluidModel', 'var') 
     mu0 = fluidModel.mu0;
     nmu = fluidModel.nmu;
     Kc  = fluidModel.Kc;
     Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
     if nmu==1, Kbc = 0; end
  else
     mu0 = 100*centi*poise;
     nmu = 0.25;
     Kc  = .1;
     Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
  end

  %% Initial vertical equilibrium
  gravity reset on, g = norm(gravity);
  [z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
  equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
  p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);
  clear equil z_0 z_max;
  
  
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
  

  %% Constant for the simulation
  dt       = totTime / numSteps;

  %% Flow equations
  phiK  = rock.perm.*rock.poro;
  gradz = grad(G.cells.centroids(:,3));
  v     = @(p, eta) ...
          -(T./(mu0*favg(eta))).*( grad(p) - g*favg(rho(p)).*gradz );
  etaEq = @(p, eta) ...
          eta - ( 1 + Kbc* cavg(v(p,eta)).^2 ./phiK ).^((nmu-1)/2);
  presEq= @(p, p0, eta, dt)  ...
          (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) + div(favg(rho(p)).*v(p, eta));


  % Define well equations
  wc = W(1).cells; % connection grid cells
  WI = W(1).WI;    % well-indices
  dz = W(1).dZ;    % connection depth relative to bottom-hole

  p_conn = @(bhp) ...
     bhp + g*dz.*rho(bhp);
  q_conn = @(p, eta, bhp) ...
     WI .* (rho(p(wc)) ./ (mu0*eta(wc))) .* (p_conn(bhp) - p(wc));
  rateEq = @(p, eta, bhp, qS) ...
     qS - sum(q_conn(p, eta, bhp))/rhoS;
  ctrlEq = @(bhp) ...
     bhp - pwf;

  %% Initialize for solution loop
  nc = G.cells.num;
  [p_ad, eta_ad, bhp_ad, qS_ad] = ...
     initVariablesADI(p_init, ones(nc,1), p_init(wc(1)), 0);
  [pIx, etaIx, bhpIx, qSIx] = ...
     deal(1:nc, nc+1:2*nc, 2*nc+1, 2*nc+2);
  sol = repmat(struct('time',[],'pressure',[],'eta',[], ...
                      'bhp',[],'qS',[]), [numSteps+1,1]);
  sol(1) = struct('time', 0, 'pressure', value(p_ad), ...
                  'eta', value(eta_ad), ...
                  'bhp', value(bhp_ad), 'qS', value(qS_ad));
  [etamin, etawmin, etamean] = deal(zeros(numSteps,1));

  %% Time loop
  t = 0; step = 0;
  while t < totTime

     % Increment time
     t = t + dt;
     step = step + 1;
     fprintf('Time step %d: Time %.2f -> %.2f days\n', ...
        step, convertTo(t - dt, day), convertTo(t, day));

     % Main Newton loop
     p0  = value(p_ad); % Previous step pressure
     [resNorm,nit] = deal(1e99, 0);
     while (resNorm > tol) && (nit < maxits)

        % Newton loop for eta (effective viscosity)
        [resNorm2,nit2] = deal(1e99, 0);
        eta_ad2 = initVariablesADI(eta_ad.val);
        while (resNorm2 > tol) && (nit2 <= maxits)
           eeq = etaEq(p_ad.val, eta_ad2);
           res = eeq.val;
           eta_ad2.val = eta_ad2.val - (eeq.jac{1} \ res);

           resNorm2 = norm(res);
           nit2     = nit2+1;
        end
        if nit2 > maxits
           error('Local Newton solves did not converge')
        else
           eta_ad.val = eta_ad2.val;
        end

        % Add source terms to homogeneous pressure equation:
        eq1     = presEq(p_ad, p0, eta_ad, dt);
        eq1(wc) = eq1(wc) - q_conn(p_ad, eta_ad, bhp_ad);

        % Collect all equations
        eqs = {eq1, etaEq(p_ad, eta_ad), ...
           rateEq(p_ad, eta_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};

        % Concatenate equations and solve for update:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
        upd = -(J \ res); % Newton update
        % Update variables
        p_ad.val   = p_ad.val   + upd(pIx);
        eta_ad.val = eta_ad.val + upd(etaIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val  = qS_ad.val  + upd(qSIx);

        resNorm = norm(res);
        nit     = nit + 1;
     end

   %  clf,
   %  plotCellData(G,eta_ad.val,'FaceAlpha',.3,'EdgeAlpha', .1);
   %  view(3); colorbar; drawnow

     if nit > maxits
        error('Newton solves did not converge')
     else % store solution
        sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
           'eta', value(eta_ad), ...
           'bhp', value(bhp_ad), 'qS', value(qS_ad));
     end
     etamin (step) = min(eta_ad.val);
     etawmin(step) = min(eta_ad.val(wc));
     etamean(step) = mean(eta_ad.val);
        
  end
endif



if switch_to=='face_based'
  %% Define geometric quantitites
  G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
  G = computeGeometry(G);

  % Discrete operators
  N  = double(G.faces.neighbors);
  intInx = all(N ~= 0, 2);
  N  = N(intInx, :);
  n = size(N,1);
  C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
  grad = @(x)C*x;
  div  = @(x)-C'*x;
  avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

  %% Rock model and transmissibilities
  pv_r = poreVolume(G, rock);
  pv   = @(p) pv_r .* exp( cr * (p - p_r) );

  hT = computeTrans(G, rock);
  cf = G.cells.faces(:,1);
  nf = G.faces.num;
  T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
  T  = T(intInx);
  fa = G.faces.areas(intInx);

  %% Basic fluid model
  rho   = @(p) rho_r .* exp( c * (p - p_r) );
  if exist('fluidModel', 'var') 
     mu0 = fluidModel.mu0;
     nmu = fluidModel.nmu;   
     Kc  = fluidModel.Kc;
     Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
     if nmu==1, Kbc = 0; end
  else
     mu0 = 100*centi*poise;
     nmu = 0.25;
     Kc  = .1;
     Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
  end

  %% Initial vertical equilibrium
  gravity reset on, g = norm(gravity);
  [z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
  equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
  p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
  

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


  %% Constant for the simulation
  dt       = totTime / numSteps;

  %% Flow equations
  phiK   = avg(rock.perm.*rock.poro).*G.faces.areas(intInx).^2;
  gradz  = grad(G.cells.centroids(:,3));
  v      = @(p, eta)   -(T./(mu0*eta)).*( grad(p) - g*avg(rho(p)).*gradz );
  etaEq  = @(p, eta)   eta - (1 + Kbc*v(p,eta).^2./phiK).^((nmu-1)/2);
  presEq = @(p, p0, eta, dt)  ...
     (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) + div(avg(rho(p)).*v(p, eta));

  %% Well model (processing specific to face-based)
  if exist('wellAvg', 'var') && wellAvg
     wavg = @(eta) 1/6*abs(C(:,W.cells))'*eta;
  else
     wavg = @(eta) ones(8,1);
  end
  
  % Define well equations
  wc = W(1).cells;
  WI = W(1).WI;
  dz = W(1).dZ;

  p_conn = @(bhp) ...
     bhp + g*dz.*rho(bhp);
  q_conn = @(p, eta, bhp) ... 
     WI .* (rho(p(wc)) ./ (mu0*wavg(eta))) .* (p_conn(bhp) - p(wc));
  rateEq = @(p, eta, bhp, qS) ...
     qS - sum(q_conn(p, eta, bhp))/rhoS;
  ctrlEq = @(bhp) ...
     bhp - pwf;

  %% Initialize for solution loop
  nc = G.cells.num;
  nf = numel(T);
  [p_ad, eta_ad, bhp_ad, qS_ad] = ...
     initVariablesADI(p_init, ones(nf,1), p_init(wc(1)), 0);
  [pIx, etaIx, bhpIx, qSIx] = ...
     deal(1:nc, nc+1:nc+nf, nc+nf+1, nc+nf+2);
  sol = repmat(struct('time',[],'pressure',[],'eta',[], ...
                      'bhp',[],'qS',[]), [numSteps+1,1]);
  sol(1) = struct('time', 0, 'pressure', value(p_ad), ...
                  'eta', value(eta_ad), ...
                  'bhp', value(bhp_ad), 'qS', value(qS_ad));
  [etamin, etawmin, etamean] = deal(zeros(numSteps,1));
  %% Time loop
  t = 0; step = 0;
  while t < totTime

     % Increment time
     t = t + dt;
     step = step + 1;
     fprintf('Time step %d: Time %.2f -> %.2f days\n', ...
        step, convertTo(t - dt, day), convertTo(t, day));

     % Main Newton loop
     p0  = value(p_ad); % Previous step pressure
     [resNorm,nit] = deal(1e99, 0);
     while (resNorm > tol) && (nit < maxits)

        % Newton loop for eta (effective viscosity)
        [resNorm2,nit2] = deal(1e99, 0);
        eta_ad2 = initVariablesADI(eta_ad.val);
        while (resNorm2 > tol) && (nit2 <= maxits)
           eeq = etaEq(p_ad.val, eta_ad2);
           res = eeq.val;
           eta_ad2.val = eta_ad2.val - (eeq.jac{1} \ res);

           resNorm2 = norm(res);
           nit2     = nit2+1;
        end
        if nit2 > maxits
           error('Local Newton solves did not converge')
        else
           eta_ad.val = eta_ad2.val;
        end

        % Add source terms to homogeneous pressure equation:
        eq1     = presEq(p_ad, p0, eta_ad, dt);
        eq1(wc) = eq1(wc) - q_conn(p_ad, eta_ad, bhp_ad);

        % Collect all equations
        eqs = {eq1, etaEq(p_ad, eta_ad), ...
           rateEq(p_ad, eta_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};

        % Concatenate equations and solve for update:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
        upd = -(J \ res); % Newton update
        % Update variables
        p_ad.val   = p_ad.val   + upd(pIx);
        eta_ad.val = eta_ad.val + upd(etaIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val  = qS_ad.val  + upd(qSIx);

        resNorm = norm(res);
        nit     = nit + 1;
     end

  %   clf,
  %   plotFaces(G,intInx, eta_ad.val,'FaceAlpha',.3,'EdgeAlpha', .1);
  %   view(3); colorbar; drawnow

     if nit > maxits
        error('Newton solves did not converge')
     else % store solution
        sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
           'eta', value(eta_ad), ...
           'bhp', value(bhp_ad), 'qS', value(qS_ad));
     end
     etamin (step) = min(eta_ad.val);
     etawmin(step) = min(wavg(eta_ad.val));
     etamean(step) = mean(eta_ad.val);
  end
endif

  
# Plot production rate and pressure decrease
addpath 'C:\Octave\Octave-5.2.0\mingw64\share\octave\5.2.0\m\plot\draw'

figure(2)
[ha,hr,hp] = ...
   plotyy([sol(2:end).time]/day, [sol(2:end).qS]*day, ...
          [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), ...
          'stairs', 'plot');
%set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
p=get(gca,'Position'); p(1)=p(1)-.01; p(2)=p(2)+.02; set(gca,'Position',p);
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