% Cut figure margin
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
perpos = get(gcf,'PaperPosition');
set(gcf,'PaperPositionMode','auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])