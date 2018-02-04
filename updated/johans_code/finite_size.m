clearvars;
hold on

flist = 1:1:8;
Ly = 50;
g = 2.0e-5;
y = 0:1:Ly;
eta_a = 1.3055;
v_th = (g*10/(2*eta_a))*y.*(y(end) - y);
grad_a = 10*g*Ly/(2*eta_a);
plot(y,v_th)
hold on
k=1;
vx = zeros(length(y), length(flist));
for i=1:length(flist)
    fname= ['Poise_vx_vy_',num2str(flist(i)),'.dat'];
    X = importdata(fname);    
    vx(:,k) = X(:,3);
    p = polyfit(y(4:end-4),vx(4:end-4,k)',2);
    eta(k) = -10*g/(2*p(1));
    grad(k) = 10*g*Ly/(2*eta(k));
    err(k) = abs(eta_a-eta(k))*100/eta_a;
    k = k+1;
end
gt = mean(vx,2);
plot(y,gt,'o')
f1 = mean(eta);
f2 = std(eta);
err_perc = abs(eta_a-f1)*100/eta_a
