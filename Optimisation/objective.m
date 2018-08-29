%EVALUATE EGGHOLDER OBJECTIVE FUNCTION

function f = objective(x, Dimension)
f=0;
for i=1:(Dimension-1)
    f = f-(x(i+1)+47)*sin(sqrt(abs(x(i+1)+0.5*x(i)+47)))-x(i)*sin(sqrt(abs(x(i)-x(i+1)-47)));
end