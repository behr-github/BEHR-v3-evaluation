function [  ] = spcaxis( clim )
%SPCAXIS( CLIM ) Changes the color limits of all subplots of a current figure

if ~isnumeric(clim) || ~isvector(clim) || numel(clim) ~= 2 || clim(2) < clim(1)
    error('spcaxis:bad_input','CLIM must be a two-element numeric vector with the first element less than the second')
end

fch = get(gcf,'children');
for a=1:numel(fch)
    if isprop(fch(a),'CLim')
        set(fch(a),'clim',clim);
    end
end


end

