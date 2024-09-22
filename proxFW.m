function [x, iter, history] = proxFW(A, b, p, sigma, opts)
% % This solves
% % min ||x||_1
% % s.t. ||Ax-b||_p <= sigma, where 1<p<2

if isfield(opts, 'L_A')
    L_A = opts.L_A;
else
    if size(A,1) < 2000
        L_A = max(eig(A*A'));
    else
        L_A = eigs(A*A',1);
    end
end

if isfield(opts, 'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = inf;
end

if isfield(opts, 'scale_beta')
    scale_beta = opts.scale_beta;
else
    scale_beta = 1;
end

if isfield(opts, 'tol')
    tol = opts.tol;
else
    tol = 5e-2;
end
if isfield(opts, 'xfeas')
    xfeas = opts.xfeas;
    x = 0*xfeas;
else
    [Q,R] = qr(A',0);
    xfeas = Q*(R'\b);
    x = 0*xfeas;
end

fval_rec = zeros(1,1);
lc_res_rec = zeros(1,1);
gap_rec = zeros(1,1);
feas_rec = zeros(1,1);

q = p/(p-1);

bd_xnm = 1 + norm(xfeas,1);

Ax = A*x;
Axb = Ax - b;
y = 0*Axb;
Axby = Axb - y;

H_0 = 1e-4;


iter = 0;
while 1
    beta_t = sqrt(iter + 1) * scale_beta;
    coefficient = beta_t / (H_0 + L_A*beta_t);
    
    % % update x
    xold = x;
    ATAxby = A'*Axby;
    xhat = x - coefficient* ATAxby;
    xplus = max(abs(xhat)-1/(H_0+ beta_t*L_A), 0) .* sign(xhat);
    x = max(min(xplus, bd_xnm), -bd_xnm);
    
    % Estimating the dual
    lc_res = norm(Axby);
    lc_res_rec(iter+1,1) = lc_res;
    lambda = beta_t*Axby;
    ATlambda = beta_t*ATAxby;
    normATlambda = norm(ATlambda,'inf');
    if normATlambda > 1
        lambda = lambda/normATlambda;
    end
    dval = -b'*lambda - sigma*norm(lambda,q);
    %%% End of computation of dual value
    
    fval = norm(x,1);
    fval_rec(iter+1,1) = fval;
    
    Ax = A*x;
    Axb = Ax - b;
    nmp_Axb = norm(Axb, p);
    feas_viol =  max(nmp_Axb - sigma, 0) ;
    feas_rec(iter+1,1) = feas_viol;
    
    % Termination criterion 1
    gap_abs = abs(fval - dval);
    gap_rel = gap_abs/max([abs(fval),abs(dval),1]);
    gap_rec(iter+1,1) = gap_rel;
    if gap_rel <= tol && feas_viol <= tol/10*sigma
        fprintf('Termination: gap and feas tolerances are satisfied: iter = %7d, fval = %6.4e, dval = %6.4e, feas = %2.1e\n', iter, fval, dval, feas_viol);
        iter = iter + 1;          % iter here is not the iteration number. This line is only used to make the length of the recorded history datas equal to iter since iter starts from 0 but the index in Matlab starts from 1.
        break;
    end
    
    
    % % update y
    yold = y;
    z = y - Axb;
    absz = abs(z);
    if sum(absz) > 1e-8
        zqp = ( sum(absz.^q) )^(1/p);
        u = ( - sigma/zqp ) * ( sign(z).* absz.^(q-1) );  % the solution of the linear oracle
        alpha = 2/(iter + 2);
        y = y + alpha * (u-y);
    else
        y = (iter/(iter+2))*y;
    end
    Axby = Axb - y;
     
    if opts.display && mod(iter, opts.verbose_freq) == 0
        fprintf('iter %7d: fval = %6.4e, lc_res = %6.2e, feas_viol = %6.2e, dval = %6.2e\n', iter, fval, lc_res,feas_viol, dval)
    end
    
    % termination criterion 2,3
    if max(norm(x-xold), norm(y-yold)) <= 1e-6
        fprintf('Termination: the linear constraint is satisfied. \n');
        iter = iter + 1;
        break;
    end
    
    if iter > maxiter
        fprintf('Termination: maximal iteration number used.\n');
        iter = iter + 1;
        break;
    end
    
    iter = iter + 1;
end
history.gap_rel = gap_rel;
history.fval_rec = fval_rec;
history.lc_res_rec = lc_res_rec;
history.gap_rec = gap_rec;
history.feas_rec = feas_rec;
