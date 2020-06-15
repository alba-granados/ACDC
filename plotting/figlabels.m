function figlabels(xlabel,ylabel,zlabel, figtitle,size,handle)
if(nargin<6)
    handle=gca;
end
    xlab = get(handle,'XLabel'); set(xlab,'String',xlabel,'FontSize',size,'FontName','Arial')
    ylab = get(handle,'YLabel'); set(ylab,'String',ylabel,'FontSize',size,'FontName','Arial')
    zlab = get(handle,'ZLabel'); set(zlab,'String',zlabel,'FontSize',size,'FontName','Arial')
    title(handle,figtitle,'FontSize',size,'FontWeight','Bold')


end