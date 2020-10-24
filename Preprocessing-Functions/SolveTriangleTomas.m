function u_next = SolveTriangleTomas(u,g,dt,h)
%Solve (I-dt*A(g))*u_next = u,where u is just a vector with boundary
%condition on two ends

sz = size(u);
d = u;
g = g*dt;
%%%%%%%%%%%slove the linear system using Tomas%%%%%%%%%%%%%%%%%%%%%%
alpha2 = ((g(1)+g(2))/(2*h^2))/(1+(g(1)+g(2))/(2*h^2));
alphan = ((g(end)+g(end-1))/(2*h^2))/(1+(g(end)+g(end-1))/(2*h^2));;
beta2 = d(1)/(1+(g(1)+g(2))/(2*h^2));
betan = d(end)/(1+(g(end)+g(end-1))/(2*h^2));

% %%%%%%%%%%%slove the linear system using Tomas%%%%%%%%%%%%%%%%%%%%%%
% alpha2 = ((g(1)+g(2))/(h^2))/(1+(g(1)+g(2))/(h^2));
% alphan = ((g(end)+g(end-1))/(h^2))/(1+(g(end)+g(end-1))/(h^2));;
% beta2 = d(1)/(1+(g(1)+g(2))/(h^2));
% betan = d(end)/(1+(g(end)+g(end-1))/(h^2));


n = length(d);
%solving the p and q
p = [alpha2];
q = [beta2];
for i=2:n-1
    a = (g(i-1)+g(i))/(2*h^2);
    c = (g(i)+g(i+1))/(2*h^2);
    b = -(a+c);
    a = -a;
    c = -c;
    b = 1-b;
    p_next = -c/(a*p(i-1)+b);
    q_next = (d(i)-a*q(i-1))/(a*p(i-1)+b);
    p = [p p_next];
    q = [q q_next];
end

% sloving the w
w = (alphan*q(end)+betan)/(1-alphan*p(end));

for i=n-1:-1:1
    w_prev = w(1)*p(i)+q(i);
    w = [w_prev w];
end

u_next = w;