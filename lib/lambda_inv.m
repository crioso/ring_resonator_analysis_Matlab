function [ n ] = lambda_inv( x , x0 )
% This function detects the index n for which x(n) = x0 holds

numberdatapoints = size(x,2);
n=[];
for i=1:(numberdatapoints-1)
    if (x(i)<= x0) && (x(i+1) >= x0)
        if (x0-x(i)) > (x(i+1)-x0)
            n=i+1;
            break;
        else
            n=i;
            break;
        end
    end
end

if isempty(n) == 1
    n=NaN;
end


