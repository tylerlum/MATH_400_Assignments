
function pde22b

' Solves the PDE using an in-built function PDEPE'

r = linspace(0,1,200);    % spatial grid for r
z = linspace(0,20,200);   % spatial grid for z
t = linspace(0,.4,161);  % output times

ra = linspace(0,1,40);    % another spatial grid
ntrunc=5;

subplot('position',[0.1 0.79 0.82 0.2])

% Calculate a Bessel function, find some zeros, and plot it all
nu=sqrt(4+1);
zng=([1:ntrunc]*2-1)*pi/2+pi/2*nu+pi/4;            % guess for zeros
for n=1:ntrunc
	zn(n)=fzero(@(z) besselj(nu,z),zng(n));    % compute zeros
end
[zn;zng]
plot(z,besselj(nu,z)),hold on
hold on,
  plot(zn,besselj(nu,zn),'*')
hold off
axis tight
xlabel('z')
ylabel('$J_{\sqrt{5}}(z)$','interpreter','latex')


%  Compute some integrals involving a Bessel function:
J=besselj(nu,sqrt(r)*zn(1));
trapz(r,2*sqrt(1-r))
trapz(r,J.^2)

% Solve PDE
sol = pdepe(0,@pdex1pde,@pdex1ic,@pdex1bc,r,t);
u = sol(:,:,1);
disp(size(u))

% Plot results
subplot('position',[0.12 0.35 0.8 0.33])
contourf(r,t,u,20)
colorbar

title('Numerical solution computed with 200 mesh points');
xlabel('Position r');
ylabel('Time t');

ns=[1 13 29 57];
t(ns)
subplot('position',[0.1 0.08 0.39 0.18])
plot(r,u(ns,:),'b');
xlabel('Position x');
ylabel('u(x,t)');
axis([0 1 -inf inf])
text(.2,.12,'t=0, 0.03, 0.07, 0.14')

subplot('position',[0.57 0.08 0.39 0.18])
plot(t,u(:,50),t,u(:,150));
xlabel('Time, t');
ylabel('u(1/4,t), u(3/4,t)');
axis([0 t(end) -inf inf])

% --------------------------------------------------------------------------

function [c,f,s] = pdex1pde(r,t,u,DuDr)

c = 1;
f = r*DuDr;
s = -u/r-DuDr;

% --------------------------------------------------------------------------

function u0 = pdex1ic(r)
  u0 = 2*sqrt(r.*(1-r));
    
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)

pl = ul;
ql = 0;
pr = ur;
qr = 0;
  