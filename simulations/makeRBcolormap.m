function cm = makeRBcolormap(varargin)

if length(varargin)==0
    lowfrac = 0.5;
else
    lowfrac=1-varargin{1};
end
A=lowfrac*100;
dA=(100-A)/63;

cm = [[ones(1,64)*100,[100:-dA:A]]',[A:dA:100,100:-dA:A]',[A:dA:100,100*ones(1,64)]'];
cm = cm/100;
cm = cm(end:-1:1,:);