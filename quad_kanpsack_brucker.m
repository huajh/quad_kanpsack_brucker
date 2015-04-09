function [ x ] = quad_kanpsack_brucker( d,a,b,b0,l,u )
%QUAD_KANPSACK_BRUCKER Summary of this function goes here

%
% x = quad_kanpsack_brucker(v, b) returns the vector x which is the solution
%   to the following constrained minimization problem:
%
%    min   sum_i^n 1/2*d_i*x_i^2 - a_i x_i
%    s.t.  sum_i^n b_i*x_i = b_0,
%           l_i <= x_i <= u_i, i=1,...,n.
%
% Author: Junhao Hua (huajh7@gmail.com)
%
%   d = [d_1,d_2,...,d_n]';
%   x = [x_1,x_2,...,x_n]';
%   a = [a_1,a_2,...,a_n]';
%   b = [b_1,b_2,...,b_n]';
%   b0 = b_0;
%   l = [l_1,l_2,...,l_n]';
%   u = [u_1,u_2,...,u_n]';
%
% Example:
%
% sigma = [0.3,0.5,0.4];
% n = 3;
% d = 2*ones(n,1);
% a = 2*sigma';
% b = 1*ones(n,1);
% l = zeros(n,1);
% u = 1*ones(n,1);
% b0 = 2;
% [ x ] = quad_kanpsack_brucker( d,a,b,b0,l,u )
%
% Reference:
%   Brucker, Peter. "An O (n) algorithm for quadratic knapsack problems."
%   Operations Research Letters 3.3 (1984): 163-166.


% critical parameters: lower and upper endpoints, tu <= tl

% dimension 
n = size(d,1);

all_tl = (a-l.*d)./b;
all_tu = (a-u.*d)./b;

% distinct critical parameter values ,Ascending
all_t = union(all_tl,all_tu);
r = size(all_t,1);

x = sol(d,a,b,l,u,all_t(1));
max_z = b'*x;
if max_z == b0
    % optimal solution
    return;
end
x = sol(d,a,b,l,u,all_t(r));
min_z = b'*x;
if min_z == b0
    % optimal solution
    return;
end

if max_z < b0 || min_z > b0
    fprintf('No feasible solution exists');
    return;
end
tmin = all_t(1);
tmax = all_t(r);
I = (1:n)'; % index

opt_num = -b0;
opt_den = 0;
while (~isempty(I))
    tl = median(all_tl(I));    
    tu = median(all_tu(I(all_tl(I)>=tl)));
    for t = [tl,tu]
        if (tmin < t && t < tmax)
           z = b'*sol(d,a,b,l,u,t);
           if z == b0
               % optimal solution     
               x = sol(d,a,b,l,u,t);
               return;       
           elseif z > b0
               tmin = max(tmin,t);
           else
               tmax = min(tmax,t);
           end      
        end
    end
    k = size(I,1);
    cache = I;
    for j=1:k
        i = cache(j);
        if all_tl(i) <= tmin
            I = setdiff(I,i);
            %x(i) = l(i);
            opt_num = opt_num + b(i)*l(i);
        end
        if tmax <= all_tu(i)
            I = setdiff(I,i);
            %x(i) = u(i);
            opt_num = opt_num + b(i)*u(i);
        end
        if all_tu(i) <= tmin && tmin <= tmax && tmax <= all_tl(i)
            I = setdiff(I,i);
            %x(i) = (a(i)-b(i)*(t/2)/d(i);
            opt_num = opt_num + b(i)*a(i)/d(i);
            opt_den = opt_den + b(i)^2/d(i);
        end        
    end
end

x = sol(d,a,b,l,u,opt_num/opt_den);
end

function x = sol(d,a,b,l,u,t)        
    % unique solution x(t)
    tmp = (a - b.*t)./d;
    x = tmp;
    id = find(tmp<=l);
    x(id) = l(id);
    id = find(tmp>=u);
    x(id) = u(id);        
end

