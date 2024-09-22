 
clear all;
rng(2)

nf = 0.01;
index_list = [4, 8, 12];
beta_list = [0.5, 1, 10, 20, 50];
repts = 20;   % repeat numbers for every beta

a = clock;

p = 1.5;    % 1<p<2
opts.display = 0;
opts.verbose_freq = 1000;
opts.maxiter = 10000;
opts.tol = 5e-2; % should only use natural choices, i.e., multiples of 5
data_all = {};

for id = 1:length(index_list)
    index = index_list(id);  % problem size
    
    %% Problem generation
    n = 2560*index;    % signal length
    k = 80*index;        % number of spikes to put down
    m = 720*index;     % number of observations to make
    
    data_betas = cell(length(beta_list) , 1);
    for repeat = 1:repts
        fprintf('            The %d-th instance\n', repeat)
        I = randperm(n);
        J = I(1:k);
        xorig = zeros(n,1);
        xorig(J) = randn(k,1);
        
        A = randn(m,n);
        for j = 1:n
            A(:,j) = A(:,j)/norm(A(:,j));
        end
        error = nf * genGauss(m, 0, 1, p);  % Generalized Gaussian noise
        b = A*xorig + error;
        sigma = 1.1*norm(error,p);
        
        [Q,R] = qr(A',0);
        xfeas = Q*(R'\b);
        opts.xfeas = xfeas;
        
        if m < 2000
            L_A = max(eig(A*A'));
        else
            L_A = eigs(A*A',1);
        end
        fprintf('End of problem generation\n')
        opts.L_A = L_A;
        
        for ib = 1:length(beta_list)
            scale_beta = beta_list(ib);
            opts.scale_beta = scale_beta;
            fprintf('scale_beta = %6.2e : \n', scale_beta)
            
            [xsoln, iter, history] = proxFW(A, b, p, sigma, opts);
            
            nnz_x = nnz(abs(xsoln)>1e-6);
            rec_err = norm(xorig - xsoln)/max(norm(xorig), 1);
            fprintf(' Recovery error %2.1e, nnz %7d\n', rec_err, nnz_x)
            
            gap_rel = history.gap_rel;
            feas_viol = history.feas_rec(end);
            data_betas{ib} = [data_betas{ib}; gap_rel, feas_viol, rec_err, nnz_x];
        end
    end
    
    matname = ['Testbeta-table-index-' int2str(index) '-' date  '-'  int2str(a(4))  '-'  int2str(a(5))  '.mat'];
    data_all{id} = data_betas;
    save(matname,'beta_list', 'data_all');
end

% % boxplots
data_index4 = data_all{1};
data_index8 = data_all{2};
data_index16 = data_all{3};


colormap = [1 0 1;
                    0 1 0;
                    0 0 1;
                    1 1 0;
                    0 1 1];

gap_beta1 = [data_index4{1}(:, 1), data_index8{1}(:, 1), data_index16{1}(:, 1) ];
gap_beta2 = [data_index4{2}(:, 1), data_index8{2}(:, 1), data_index16{2}(:, 1) ];
gap_beta3 = [data_index4{3}(:, 1), data_index8{3}(:, 1), data_index16{3}(:, 1) ];
gap_beta4 = [data_index4{4}(:, 1), data_index8{4}(:, 1), data_index16{4}(:, 1) ];
gap_beta5 = [data_index4{5}(:, 1), data_index8{5}(:, 1), data_index16{5}(:, 1) ];

% % plot the gap for all beta's
gap_plot = {gap_beta1; gap_beta2; gap_beta3;gap_beta4; gap_beta5};
fig1 = figure;
labels = ['i= 4'; 'i= 8';'i=12'];
aboxplot(gap_plot,'labels', labels, 'colormap', colormap);
ylabel('${\rm gap}_r$','Interpreter','latex')
legend('$\beta_0$=0.5','$\beta_0$=1','$\beta_0$=10', '$\beta_0$=20', '$\beta_0$=50','Interpreter', 'latex',...
    'Location','northwest')
% title('gap with different beta')
a = clock;
figname = ['Gap-boxplot-date-' num2str(a(4)) '-'  num2str(a(5))];
savefig(fig1, figname, 'Figures');


% % plot the feas for all beta's
feas_beta1 = [data_index4{1}(:, 2), data_index8{1}(:, 2), data_index16{1}(:, 2) ];
feas_beta2 = [data_index4{2}(:, 2), data_index8{2}(:, 2), data_index16{2}(:, 2) ];
feas_beta3 = [data_index4{3}(:, 2), data_index8{3}(:, 2), data_index16{3}(:, 2) ];
feas_beta4 = [data_index4{4}(:, 2), data_index8{4}(:, 2), data_index16{4}(:, 2) ];
feas_beta5 = [data_index4{5}(:, 2), data_index8{5}(:, 2), data_index16{5}(:, 2) ];

feas_plot = {feas_beta1; feas_beta2; feas_beta3;feas_beta4;feas_beta5};

fig2 = figure;
labels = ['i= 4'; 'i= 8';'i=12'];
aboxplot(feas_plot,'labels', labels, 'colormap', colormap); % Advanced box plot
%
% xlabel(''); % Set the X-axis label
ylabel('Feasiblity violations')
legend('$\beta_0$=0.5','$\beta_0$=1','$\beta_0$=10', '$\beta_0$=20', '$\beta_0$=50','Interpreter','latex',...
    'Location','northwest')
% title('feasiblility violation with different beta')

figname = ['Feasibility-boxplot-date-' num2str(a(4)) '-'  num2str(a(5))];
savefig(fig2, figname, 'Figures');
