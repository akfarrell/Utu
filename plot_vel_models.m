% ------ Load Data -------- %
vel_s = load('Heather_vel_things/3D_1.75_Utur_0.1km.mat') %structure of velocities
vel_deps = vel_s.zmod;
p_vels = vel_s.apvel;
vel_deps_corr = (-vel_deps+3.5); 
p_vels_extrap = p_vels;
p_vels_extrap(76:329) = p_vels(75); %use 5712.2447 for values that would be in magma body


% ------ Plot Figure -------- %
close all;
r = figure;
set(r, 'Position', [1000 1000 1000 1000])
plot(p_vels,vel_deps_corr,'k--')
hold on
plot(p_vels_extrap, vel_deps_corr, 'k')
%set(gca,'yticklabel',num2str(get(gca,'ytick')')) %to keep y in
%non-scientific notation
xlabel({'';'P-wave velocity (km/s)'}, 'FontSize', 14)
ylabel('Depth (km)', 'FontSize', 14)
ylim([vel_deps_corr(end)-1 vel_deps_corr(1)+1]);
set(gca,'fontsize', 14)
hgexport(r, 'Velocities_used.png', hgexport('factorystyle'), 'Format', 'png');