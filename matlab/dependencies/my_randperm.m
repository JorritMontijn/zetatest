function x = my_randperm(n,k)
% randperm introduced to make results reproducable between python and
% MATLAB implementation

if nargin < 2
    k = n;
end

[~,ind] = sort(rand(n,1));
x = ind(1:k);
end