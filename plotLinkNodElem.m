function plotLinkNodElem(nod, elem)
[num_nod, dim]=size(nod);
num_elem=size(elem,1);
xmax=max(nod(:,1));
xmin=min(nod(:,1));
if (dim == 1)
    ymax=0.5;
    ymin=-0.5;
else 
    ymax=max(nod(:,2));
    ymin=min(nod(:,2));
end
shift=max(0.02*max(abs([xmax, xmin,ymin,ymax])),0.01);
if (dim == 1) %Plot 1D nodes
    plot(nod(:,1),0,'ro','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10); 
    axis([xmin-shift, xmax+shift, ymin-shift, ymax+shift])
    hold on;
    for j=1:num_nod   
        text(nod(j,1), -0.5*shift,['n' num2str(j,'%2d')]);
    end
    plot(nod(:,1),zeros(num_nod,1),'-','LineWidth',2); %Plot elements
    hold on;
    for e=1:num_elem   
        text(0.5*(nod(e,1)+nod(e+1,1)), shift,['E' num2str(e,'%2d')]);
    end
else
    plot(nod(:,1),nod(:,2),'ro','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
    axis([xmin-shift, xmax+shift, ymin-shift, ymax+shift])
    hold on;
    for j=1:num_nod   
        text(nod(j,1), nod(j,2)-shift,['n' num2str(j,'%2d')]);
    end
    %Plot elements
    for e=1:num_elem  
        plot([nod(elem(e,1),1); nod(elem(e,2),1)],...
             [nod(elem(e,1),2); nod(elem(e,2),2)],'-b','LineWidth',2);
        text(0.5*(nod(elem(e,1),1)+nod(elem(e,2),1))-shift, ...
             0.5*(nod(elem(e,1),2)+nod(elem(e,2),2))-shift,['E' num2str(e,'%2d')]);
    end
end
hold off;

