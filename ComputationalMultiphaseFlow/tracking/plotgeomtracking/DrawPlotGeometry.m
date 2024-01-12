function DrawPlotGeometry(yLim, R_geom)
%% Draw Plot Geometry
    bigcircle = [-.5 -.5 1 1]; bigrectangle = [0 -0.5 1 1];
    figure
    hold on
    subplot(1,2,1)
    rectangle('Position', bigrectangle,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9],'LineWidth',2);
    xline1 = [0 1]; yline1 = [0 0];
    line(xline1,yline1, 'LineWidth',2, 'Color', [0 1 0], 'LineStyle','--');
    xline2 = [0 1]; yline2 = [0.5 0.5];
    line(xline2,yline2,'LineWidth',2,'Color',[1 0 0]);
    xline3 = [0 0]; yline3 = [yLim/(2*R_geom) 0.5];
    line(xline3,yline3,'LineWidth',2,'Color',[1 0 0]);
    xline4 = [0 1]; yline4 = [-0.5 -0.5];
    line(xline4,yline4,'LineWidth',2,'Color',[1 0 0]);
    xline5 = [0 0]; yline5 = [-yLim/(2*R_geom) -0.5];
    line(xline5,yline5,'LineWidth',2,'Color',[1 0 0]);
    hold on
    subplot(1,2,2)
    rectangle('Position', bigcircle,'Curvature',[1 1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','red','LineWidth',2);
end