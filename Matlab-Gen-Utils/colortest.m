function [  ] = colortest( colorin )
% colortest(colorspec): Shows what a color will look like

figure; fill([1 1 2 2],[1 2 2 1],colorin);
if ~ischar(colorin);
    s = mat2str(colorin);
else
    s = colorin;
end
    title(s,'fontsize',18);

end

