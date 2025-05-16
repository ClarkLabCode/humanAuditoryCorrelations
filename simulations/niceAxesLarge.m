%% makes current axes nice

set(gca,'fontname','times','fontsize',20);
ax=findall(gca,'Type','line');
for i=1:length(ax)
    set(ax(i),'linewidth',2);
end
ax=findall(gcf,'Type','text');
for i=1:length(ax)
    set(ax(i),'fontname','times','fontsize',20);
end