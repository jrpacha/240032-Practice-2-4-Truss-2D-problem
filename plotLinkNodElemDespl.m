function plotLinkNodElemDespl(nod, elem, u, esc)
[num_nod, dim]=size(nod);
num_elem=size(elem,1);
nodu=zeros(num_nod,dim);
if (dim == 1)
    ymax=0.5;
    ymin=-0.5;
    nodu=nod+u;
else 
    nodu(:,1)=nod(:,1)+esc*u(1:2:num_nod*dim,1);
    nodu(:,2)=nod(:,2)+esc*u(2:2:num_nod*dim,1);
    ymax=max(nodu(:,2));
    ymin=min(nodu(:,2));

end
xmax=max(nodu(:,1));
xmin=min(nodu(:,1));

shift=max(0.02*max(abs([xmax, xmin,ymin,ymax])),0.01);
if (dim == 1) %Plot nodes
    plot(nodu(:,1),0,'ro','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); 
    %axis([xmin-shift, xmax+shift, ymin-shift, ymax+shift])
    hold on;
    plot(nod(:,1),zeros(num_nod,1),'-','LineWidth',1); %Plot elements
    plot(nodu(:,1),zeros(num_nod,1),'-','LineWidth',1); %Plot elements desp
else
    axis([xmin-shift, xmax+shift, ymin-shift, ymax+shift])
    hold on;
    %Plot elements
    for e=1:num_elem  
        plot([nod(elem(e,1),1); nod(elem(e,2),1)],...
             [nod(elem(e,1),2); nod(elem(e,2),2)],'b-.','LineWidth',1);
        plot([nodu(elem(e,1),1); nodu(elem(e,2),1)],...
             [nodu(elem(e,1),2); nodu(elem(e,2),2)],'-b','LineWidth',1);
    end
end
hold off;

