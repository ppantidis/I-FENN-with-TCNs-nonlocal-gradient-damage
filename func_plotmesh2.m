function func_plotmesh2(coords,connect,nelem,damage,elocal,enonlocal,nelnodes,color,model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================== PLOTTING ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust nonlocal strain graph c-axis limits
% caxis_minval = min(min(elocal), min(enonlocal));
% caxis_maxval = max(max(elocal), max(enonlocal));

% Initialize figure
f = figure("Position",[400 200 750 400]);
t = tiledlayout(2,2,'TileSpacing','Compact');
title(t, strcat(model_name + ": Increment #" + num2str(increment) + " - Iteration #" + num2str(last_iteration) + " - Loadfactor " + num2str(loadfactor)), 'fontweight','bold');


% -------------------------------------------------------------------------
% ------------------------ PLOTTING DAMAGE CONTOUR ------------------------
% -------------------------------------------------------------------------
nexttile
hold on; box on;
colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
f2D_4 = [1,2,3,4];
caxis([0 1])
for lmn = 1:nelem
   for i = 1:nelnodes(lmn)
       x(i,1:2) = coords(1:2,connect(i,lmn));
       col(i,1) = damage(connect(i,lmn));
   end

   if (nelnodes(lmn)==4)
       patch('Vertices',x,'Faces',f2D_4,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','k');
   else
       disp("Check your nelnodes")
   end
end
colorbar
title("Damage contours", 'fontweight','bold');
axis equal
hold off

% -------------------------------------------------------------------------
% -------------------------- PLOTTING RESIDUALS ---------------------------
% -------------------------------------------------------------------------
nexttile
semilogy(0:last_iteration,Res_F_F_norm,":ro",'LineWidth',1.1);
hold on; box on; grid on;
semilogy(1:last_iteration,Res_u_norm,"--bs",'MarkerFaceColor','b');
xlim([0, last_iteration]); 
set(gca, 'XTick', 1:last_iteration); xlabel("No. iterations"); 
legend(["Internal stresses Residual", "Displacements Residual"]);
title("Residuals (norm)", 'fontweight','bold');
hold off



% -------------------------------------------------------------------------
% ----------------------- PLOTTING ELOCAL CONTOUR -------------------------
% -------------------------------------------------------------------------
nexttile
hold on; box on;
colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
f2D_4 = [1,2,3,4];
%caxis([0 5e-3]);
for lmn = 1:nelem
   for i = 1:nelnodes(lmn)
       x(i,1:2) = coords(1:2,connect(i,lmn));
       col(i,1) = elocal(connect(i,lmn));
   end

   if (nelnodes(lmn)==4)
       patch('Vertices',x,'Faces',f2D_4,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','none');
   else
       disp("Check your nelnodes")
   end
end
colorbar
title("Local strain contours", 'fontweight','bold');
axis equal
hold off


% -------------------------------------------------------------------------
% --------------------- PLOTTING ENONLOCAL CONTOUR ------------------------
% -------------------------------------------------------------------------
nexttile
hold on; box on;
colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
f2D_4 = [1,2,3,4];
%caxis([0 2.5e-3]);
for lmn = 1:nelem
   for i = 1:nelnodes(lmn)
       x(i,1:2) = coords(1:2,connect(i,lmn));
       col(i,1) = enonlocal(connect(i,lmn));
   end

   if (nelnodes(lmn)==4)
       patch('Vertices',x,'Faces',f2D_4,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','k');
   else
       disp("Check your nelnodes")
   end
end
colorbar
title("Nonlocal strain contours", 'fontweight','bold');
axis equal
hold off


% -------------------------------------------------------------------------
saveas(f, strcat(model_name + "_Damagecontour_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".png"))

close all 


end













