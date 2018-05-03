function Px=BoxProjection(x,bounds)
% Projection of x onto the closed convex set C={x: bounds(1)<= x <= bounds(2)}

lb=bounds(1);%lower box bound
ub=bounds(2);%upper box bound

if isequal(lb,-Inf) && isequal(ub,Inf)
  Px=x;
elseif isequal(lb,-Inf) && isfinite(ub)
  x(x>ub)=ub;
  Px=x;
elseif isequal(ub,Inf) && isfinite(lb)
  x(x<lb)=lb;
  Px=x;
else
  x(x<lb)=lb;
  x(x>ub)=ub;
  Px=x;
end