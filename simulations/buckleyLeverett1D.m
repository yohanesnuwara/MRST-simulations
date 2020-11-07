###########################################################################
#    Two-phase Incompressible and Immiscible Simulation in 1D Reservoir   #
#        Waterflooding in Dead (non-volatile; non-gas-saturated) Oil      #
#                                                                         #
#                      Buckley-Leverett Simulation                        #
###########################################################################

% We investigate the effect of time step on the accuracy and convergence of
% the implicit transport solver. We find that 200 timesteps for implicit is
# the best. See explanations in p. 299

mrstModule add incomp

G = computeGeometry(cartGrid([100,1]));
rock = makeRock(G, 100*milli*darcy, 0.2);

fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000, 1000] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
bc = fluxside([], G, 'Left',   1, 'sat', [1 0]);
bc = fluxside(bc, G, 'Right', -1, 'sat', [0 1]);

hT = computeTrans(G, rock);
rSol = initState(G, [], 0, [0 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'bc', bc);

% Explicit tranport solver
%rSole = explicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'verbose', true);

% Implicit transport solver: try with one time step
[rSoli, report] = ...
    implicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'Verbose', true);

%% Plot results with various number of time steps
figure(1)
plot(G.cells.centroids(:,1), rSole.s(:,1),'k--','LineWidth',1.5);
leg = cell(7,1); leg{1} = 'Expl: 199 steps';

n   = [4 10 20 40 100 200];
its = [0 0 0 0 0 0 0];
col = 'rgbcmk';
hold on
for k=1:numel(n)
    rSolt = rSol;
    for i=1:n(k)
        [rSolt, report] = ...
            implicitTransport(rSolt, G, 10/n(k), rock, fluid, 'bc', bc);
        its(k) = its(k) + report.iterations + report.vasted_iterations;
    end
    plot(G.cells.centroids(:,1),rSolt.s(:,1), [col(k) '-'],'LineWidth',1.5);
    leg{k+1} = sprintf('n=%3d: %3d its',n(k),its(k));
end
hold off
legend(leg{:});
title('Saturation Distribution for Various Number of Timesteps')

%% Plot the best result for maximum number of timesteps = 200
figure(2)
plot(G.cells.centroids(:,1),rSolt.s(:,1), [col(k) '-'],'LineWidth',1.5);
title('Saturation Distribution: Best Result using n=200')