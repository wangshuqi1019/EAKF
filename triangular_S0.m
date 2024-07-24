% version from Zengyang
function [R] = triangular_S0(m,l,u)
%TRIANGULAR_S0 Summary of this function goes here
%   generate 100 S0 by triangular distribution with mean, lower_95,
%   upper_95
%     syms a b c
%     e1 = (l-a)^2/((c-a)*(b-a));
%     e2 = (b-u)^2/((b-c)*(b-a));
%     e3 = (a+b+c)/3;
%     eqs = [e1==0.05, e2==0.05, e3==m, a<b, c<b, a<c,a<l,b>u];
%     res = solve(eqs,a,b,c);
% 
%     pd = makedist("Triangular","a",double(res.a),"b",double(res.c),"c",double(res.b));
%     R = random(pd,[100,1]);
    
    b = m;
    a=l;
    c=u;

    pd = makedist("Triangular","a",a,"b",b,"c",c);
    R = random(pd,[100,1]);
end

