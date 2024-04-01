function func_plotmesh_gps(xcoord,ycoord,damage,elocal,enonlocal_pred,model_name,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================== PLOTTING ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust nonlocal strain graph c-axis limits
% caxis_minval = min(min(elocal), min(enonlocal));
% caxis_maxval = max(max(elocal), max(enonlocal));
scattersize = 6;
geom_xlimit = [-5 105];
geom_ylimit = [-5 105];
caxismax = 1e-6 * loadfactor / 0.005;

% Initialize figure
f = figure("Position",[400 200 650 550]);
t = tiledlayout(2,2,'TileSpacing','Compact');
title(t, strcat(model_name + ": Increment #" + num2str(increment) + " - Iteration #" + num2str(last_iteration) + " - Loadfactor " + num2str(loadfactor)), 'fontweight','bold');


% -------------------------------------------------------------------------
% ------------------------ PLOTTING DAMAGE CONTOUR ------------------------
% -------------------------------------------------------------------------

nexttile(1); hold on; box on; set(gca,'ytick',[0 20 40 60 80 100]); set(gca,'xtick',[0 20 40 60 80 100]); 
scatter(xcoord,ycoord,scattersize,damage,'filled');
xlim(geom_xlimit); ylim(geom_ylimit); title1 = title('\boldmath${\overline{d}_{pred}}$'); colormap(gca,jet); colorbar
set(gca,'fontsize',10); set(title1,'fontsize',11); 
caxis([0 1]);
hold off

% -------------------------------------------------------------------------
% -------------------------- PLOTTING RESIDUALS ---------------------------
% -------------------------------------------------------------------------
nexttile(2); 
semilogy(0:last_iteration,Res_F_F_norm,":ro",'LineWidth',1.1);
hold on; box on; grid on;
semilogy(1:last_iteration,Res_u_norm,"--bs",'MarkerFaceColor','b');
xlim([0, last_iteration]); 
set(gca, 'XTick', 1:last_iteration); xlabel("No. iterations"); 
legend(["Int stress", "Disp"]);
title("Residuals (norm)", 'fontweight','bold');
hold off



% -------------------------------------------------------------------------
% ----------------------- PLOTTING ELOCAL CONTOUR -------------------------
% -------------------------------------------------------------------------

nexttile(3); hold on; box on; set(gca,'ytick',[0 20 40 60 80 100]); set(gca,'xtick',[0 20 40 60 80 100]); 
scatter(xcoord,ycoord,scattersize,elocal,'filled');
xlim(geom_xlimit); ylim(geom_ylimit); title1 = title('\boldmath${{\epsilon}_{pred}}$'); colormap(gca,jet); colorbar
set(gca,'fontsize',10); set(title1,'fontsize',11);
hold off


% -------------------------------------------------------------------------
% --------------------- PLOTTING ENONLOCAL CONTOUR ------------------------
% -------------------------------------------------------------------------

nexttile(4); hold on; box on; set(gca,'ytick',[0 20 40 60 80 100]); set(gca,'xtick',[0 20 40 60 80 100]); 
scatter(xcoord,ycoord,scattersize,enonlocal_pred,'filled');
xlim(geom_xlimit); ylim(geom_ylimit); title1 = title('\boldmath${\overline{\epsilon}_{pred}}$'); colormap(gca,jet); colorbar
set(gca,'fontsize',10); set(title1,'fontsize',11);
caxis([0 caxismax])
hold off



% -------------------------------------------------------------------------
saveas(f, strcat(model_name + "_Damagecontour_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".png"))

close all 


end













