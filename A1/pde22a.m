function pde22a

' Solves the PDE using an in-built function PDEPE'

m = 0;
x = linspace(0,pi/2,200);   % spatial grid
t = sort([ .02 .08 .2 .8 2 linspace(0,1,101).^2*4])     % output times

N = 1000;
[X, T] = meshgrid(x, t);
u = zeros(numel(t), numel(x));

for n = 1:N
    u(:, :) = u(:, :) + c(n) * cos((2*n - 1)*X) .* exp(-(2*n - 1)^2 .* T);
end

subplot('position',[0.1 0.48 0.88 0.46])
%imagesc(x,t,u),axis xy
%contour(x,t,u,20)
contourf(x,t,u,20)
  colorbar
%pcolor(x,t,u),shading interp

title('20-term series solution computed with 200 mesh points','fontsize',12);
xlabel('Position x','fontsize',12);
ylabel('Time t','fontsize',12);

set(gca,'yscal','log')

ns=[1 9 17 26 49 76]; t(ns)
subplot('position',[0.1 0.09 0.4 0.29])
plot(x,u(ns,:),'b');
xlabel('Position x','fontsize',12);
ylabel('u(x,t)','fontsize',12);
axis([0 pi/2 0 .63])
text(.2,.68,'t=0, 0.02, 0.08, 0.2, 0.8, 2','fontsize',12)


subplot('position',[0.58 0.09 0.4 0.29])
plot(t.^.5,2/pi*(4-pi)*exp(-t),'r--','linewidth',1);
hold on
plot(t.^.5,sqrt(2)/pi*(4-pi)*exp(-t),'r--','linewidth',1);
plot(t.^.5,u(:,1),'b')
plot(t.^.5,u(:,100),'color',[0 .5 0])
xlabel('(Time,t)^{1/2}','fontsize',12);
ylabel('u(0,t) and u(\pi/4,t)','fontsize',12);
axis([0 t(end).^.5 0 inf])
text(1,.24,'\leftarrow u(0,t)','fontsize',13,'color','b')
text(.5,.1,'u(\pi/4,t) \rightarrow','fontsize',13,'color',[ 0 .5 0])


% --------------------------------------------------------------------------

function c = c(n)
if mod(2*n - 1, 4) == 1
    c = (8 - 2*pi*(2*n - 1)) / (pi * (2*n - 1)^3);
else
    c = (-8 - 2*pi*(2*n - 1)) / (pi * (2*n - 1)^3);
end
