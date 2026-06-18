function [ x ] = make_long( x )
%just makes sure a variable is long or wide
if size(x,2) > 1
    x=x';
end
return




















