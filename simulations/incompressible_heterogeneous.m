###################################################################
# Incompressible (Water) Simulation on 3D Heterogeneous Reservoir # 
###################################################################

mrstModule add incomp

# Define geometry
nx = 10; ny = 10; nz = 10;
G = cartGrid([nx ny nz]);
%G.nodes.coords = twister(G.nodes.coords); # activate this to get twisted cells
G = computeGeometry(G);

# Distribute Gaussian on porosity
p = gaussianField(G.cartDims, [0.2 0.4], [5 3 1], 2.5);

# Distribute logNorm on NTG
ntg = logNormLayers(G.cartDims, [0.2 0.4]);

# Transform to permeability using Kozeny-Carmen
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);

# Make rock based on poro, perm, and NTG
rock = makeRock(G, K(:), p(:), 'ntg', ntg);

# Define fluid
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

# Initialize solver                            
resSol = initResSol(G, 0.0);

# Impose boundary conditions
bc = fluxside([], G, 'LEFT',  100*meter^3/day());
bc = fluxside(bc, G, 'FRONT', 50*barrel()/day());
bc = pside   (bc, G, 'RIGHT', 100*barsa);

# Give wells
%% add vertical well
%% producer 1 with constant rate 10 bbl/d
cellsWell1 =  1 : nx*ny : nx*ny*nz;
radius     = .1; % [m]
W = addWell([], G, rock, cellsWell1, 'Type', 'rate', ...
            'InnerProduct', 'ip_tpf', ...
            'Val', -10.*barrel()/day(), 'Radius', radius, 'name', 'I');

%% add horizontal well in y-direction
%% injector well 2 with constant pressure
cellsWell2 = nx : ny : nx*ny;
W = addWell(W, G, rock, cellsWell2, 'Type', 'bhp', ...
            'InnerProduct', 'ip_tpf', ...
            'Val', 110*barsa, 'Radius', radius, 'Dir', 'y', 'name', 'P');

# Compute transmissibility
T = computeTrans(G, rock, 'Verbose', true);

# Solve pressure using TPFA
resSol = incompTPFA(resSol, G, T, fluid, ...
                   'bc', bc, ...
                   'wells', W, ...
                   'MatrixOutput', true);

# Plot outputs
clf

subplot(2,2,1)
plotCellData(G, convertTo(rock.perm,milli*darcy));                     
title('Cell Permeability [md]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); shading faceted; 
axis tight;
colorbar   

subplot(2,2,2)
plotCellData(G,rock.poro,'EdgeColor','none');                     
title('Cell Porosity [v/v]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); shading faceted; 
axis tight;
colorbar   

subplot(2,2,3)
plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()), ...
             'EdgeColor', 'none');                     
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); shading faceted; 
axis tight;
colorbar            

%{
# plot wells (NOTE: not working properly)
subplot(1,2,1)
plotCellData(G,rock.poro,'EdgeColor','none');
axis equal tight; view(3);

subplot(1,2,2)
plotCellData(G,convertTo(rock.perm,milli*darcy), 'EdgeColor', 'none');
axis equal tight; view(3);
%}