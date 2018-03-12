set(gcf,'Units','Inches');
pos = get(gcf,'Position');
perpos = get(gcf,'PaperPosition');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])