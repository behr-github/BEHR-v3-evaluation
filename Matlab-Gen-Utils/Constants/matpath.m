function [ path ] = matpath(  )
%matpath Returns userpath without the trailing colon

path = userpath;
path = path(1:end-1);

end

