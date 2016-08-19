% Visualize a permeability layer of SPE10 project
% previously saved in MATLAB format 
load 'spe_layer1_perm.mat' perm_layer;
rock.perm = perm_layer*9.869232667160130e-13;

figure;
rock.k = rock.perm(:,1);
kk =  reshape(log10(rock.k),60,220)';
contourf(linspace(0,10,60), linspace(0,10,220), kk, 21 )
axis tight, colorbar
drownow;
