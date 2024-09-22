
% clear all;
rng(2)

nf = 0.01;
index = 8; % problem size

%% Problem generation
p = 1.5;
n = 2560*index;      % signal length
k = 80*index;          % number of spikes to put down 
m = 720*index;       % number of observations to make

I = randperm(n);
J = I(1:k);
xorig = zeros(n,1);
xorig(J) = randn(k,1);

A = randn(m,n);
for j = 1:n
    A(:,j) = A(:,j)/norm(A(:,j));    
end
error = nf*genGauss(m, 0, 1, p);     % Generalized gaussian noise; the last input parameter control the shape of the distribution
b = A*xorig + error;
sigma = 1.1*norm(error,p);

[Q,R] = qr(A',0);
xfeas = Q*(R'\b);
opts.xfeas = xfeas;

tstart1 = tic;
if m < 2000
    L_A = max(eig(A*A'));
else
    L_A = eigs(A*A',1);
end
time_eigs = toc(tstart1);

fprintf('End of problem generation\n')
opts.L_A = L_A;
opts.display = 1;
opts.verbose_freq = 1000;
opts.maxiter = 10000;
opts.scale_beta = 20;     
opts.tol = 5e-2;     % should only use natural choices, i.e., multiples of 5

tstart = tic;
[xsoln, iter, history] = proxFW(A, b, p, sigma, opts);
time_proxFW = toc(tstart);

nnz_x = nnz(abs(xsoln)>1e-6);
rec_err = norm(xorig - xsoln)/max(norm(xorig), 1);
fprintf(' Recovery error %2.1e, nnz %7d\n', rec_err, nnz_x)

a = clock;
matname = ['history_index-' num2str(index) '-beta0-'  num2str(opts.scale_beta)  '-' date '-' num2str(a(4)) '-'  num2str(a(5)) '.mat'];
save(matname,'xorig', 'xsoln', 'history', 'iter', 'index');

fig1 = figure;
plot(1:n, xorig, 'bo', 1:n, xsoln, 'r*')
xlim([0,n])

figname = ['recovery_plot_index_'  num2str(index) '-' date '-' num2str(a(4)) '-'  num2str(a(5))];
savefig(fig1, figname, 'Figures');
% % 

fig2 = figure;
loglog(1:iter-1, history.lc_res_rec(2:iter), 'r-', 'LineWidth', 1.5)  % % the entry lc_res_rec(iter) stores the info at x^{iter-1} since we start from iter=0 but Matlab index start from 1
hold on 
loglog(1:iter-1,history.gap_rec(2:iter), 'm-', 'LineWidth',1.5)
hold on
loglog(1:iter-1,1./sqrt(2:iter), 'b-', 'LineWidth', 1.5)      % here we plot 1/sqrt(t+1)
% legend({'$\|Ax-b-y\|$', '$|f(x^k)-f_{end}|$', '$1/{\sqrt{k}}$'}, 'Interpreter','latex')
legend({'$\|Ax^t-b-y^t\|$', '${\rm gap}_r(t)$', '$1/{\sqrt{t+1}}$'}, 'Interpreter','latex','Location', 'NorthEast')
xlabel('$t$','Interpreter','latex')

figname = ['complexity_plot_index_'  num2str(index) '-' date '-' num2str(a(4)) '-'  num2str(a(5))];
savefig(fig2, figname, 'Figures');
