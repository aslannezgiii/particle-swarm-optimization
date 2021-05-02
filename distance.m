
function z = distance(x,target)

z = sqrt(sum(((x(1,1)-target(1,1)).^2) + ((x(1,2)-target(1,2)).^2) + ((x(1,3)-target(1,3)).^2)) );

end
