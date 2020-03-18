G = cartGrid([1, 1, 30], [1, 1, 30]*meter^3);
G = computeGeometry(G);
# plotGrid(G); view(3);

## rock with 100 md perm
rock = makeRock(G, 0.1*darcy(), 0.2);
cla, plotCellData(G, rock.perm)

## a structure to hold the reservoir state
sol = initResSol(G, 0.0);

## pressure 100 bar at top of the column, no flow (v=0)
bc = pside([], G, 'TOP', 100.*barsa());

## transmissibility coefficient of proportionality Tij
hT = computeTrans(G, rock);

## fluid properties: incompressible flow
mrstModule add incomp
gravity reset on
# fluid with density 1014 kg/m3 & 1 cp, equals to water
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
[mu, rho] = fluid.properties();

# incompressible solver
sol = incompTPFA(sol, G, hT, fluid,'bc', bc);

# plot the solution
plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
set(gca, 'ZDir', 'reverse'), title('Pressure [bar]'), colormap(jet), caxis([100 103])
view(3), colorbar, set(gca,'DataAspect',[1 1 10])
