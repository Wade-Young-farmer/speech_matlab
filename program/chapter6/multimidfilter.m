function y=multimidfilter(x,m)
a=x;
% m refers to the times of doing med filtering
for k=1 : m
    b=medfilt1(a, 5); 
    a=b;
end
y=b;