function [ e ] = ste(M)
%
s = size(M,1)
e = std(M)./sqrt(s-1)

end

