classdef IH
methods(Static)
    %% 1: Infinite Horizon models:
    
    % 1: Define griddedInterpolant for labor equilibrium, as function of
    % tau_t and tau_{t+1}:

    function [hfunc,X1,X2,h] = h_EE_approximation(par,approx)
        % Grids:
        tau = base.nonlinspace(eps,1-eps,approx.n,approx.s);
        tau_p = base.nonlinspace(eps,1-eps,approx.nn,approx.ss);
        [X1,X2] = ndgrid(tau,tau_p);
        h0 = zeros(size(X1(1,:)))+0.5;
        % Solve for h on the grid:
        h = NaN(size(X1));
        for i=1:size(X1,1)
            h(i,:) = base.h_equi(par,X1(i,:),X2(i,:),par.epsilon,par.theta,h0);
            h0 = h(i,:);
        end
        hfunc = griddedInterpolant(X1,X2,h,'spline','linear');
    end

    % 2: Infinite horizon policy function: Gradient free solver,
    % pre-approximated PEE labor function.
    function [sample,dif,maxdiff,x,tau_cell,tau_m_cell,h] = GradientFree_Policy(par,np,grids,hfunc,tau0)
        % 1: tau0 is a first guess on a policy function
        % tau(tau_{t-1},nu_t). If nan, guess constant function:
        if isnan(tau0)
            tau0 = @(tau,nu) 0.5;
        end
        % 2: Data structure:
        maxdiff= NaN(np.max_iter,1);
        [sample,dif] = deal(cell(np.max_iter,1));
        [tau_m_cell,h] = deal(cell(np.n_nu,1));
        % 3: Iteration:
        for x=1:np.max_iter
            if x>1
                tau0 = @(tau,nu) taufunc(tau,nu);
            end
            % Step 1: Given policy guess, find optimal tau:
            tau_cell = IH.pol_gf_grids(par,np,grids.h,grids.nu,hfunc,tau0);
            % Step 2: Back out optimal tau_t on grid of tau_{t-1}:
            tau_m_vec = base.BOT(par,repmat(grids.h',np.n_nu,1),vertcat(tau_cell{:}),par.epsilon,par.theta);
            % Step 3: Remove NaN (nb: no need for vectorization here)
            for i=1:length(tau_cell)
                tau_m_cell{i} = tau_m_vec((i-1)*np.n_h+1:i*np.n_h);
                tau_cell{i} = tau_cell{i}(~isnan(tau_m_cell{i}));
                h{i} = grids.h(~isnan(tau_m_cell{i}))';
                tau_m_cell{i} = rmmissing(tau_m_cell{i});
            end
            % Step 4: Update policy function with interpolant:
            taufunc = @(tau,i) IH.update_taufunc(tau_cell,tau_m_cell,tau,i);
            
            % Step 5: Test convergence of function
            sample{x} = NaN(np.n_test,np.n_nu);
            dif{x} = NaN(np.n_test,np.n_nu);
            for i=1:np.n_nu
               sample{x}(:,i) = taufunc(grids.test',i);
               dif{x}(:,i) = taufunc(grids.test',i)-tau0(grids.test',i);
            end
            maxdiff(x) = max(abs(dif{x}),[],'all');
            if maxdiff(x)<np.tol_diff
                break
            end
            
        end
    end
  
    function taufunc = update_taufunc(tau_array,tau_m_array,tau,i)
        [tau_m,index] = sort(tau_m_array{i});
        func = griddedInterpolant(tau_m, tau_array{i}(index),'pchip');
        taufunc = func(tau);
    end

    % 3: Optimal policy with gradient free solver.
    function tax_cell= pol_gf_grids(par,np,h_grid,nu_grid,hfunc,taufunc)
        tax_cell = cell(length(nu_grid),1);
        options = optimset('TolX',np.tol_solver);
        
        for i=1:length(nu_grid)
            taufunc_nu = @(tau) taufunc(tau,i);
            hpolicy = @(tau) hfunc(tau,taufunc_nu(tau));
            gamma0_m = @(tau,hi) base.gamma0f(par,hi,tau,par.epsilon);
            gamma1_m = @(tau,hi) base.gamma1f(par,hi,tau,par.epsilon);
            gamma0 = @(tau) base.gamma0f(par, hpolicy(tau), taufunc_nu(tau), par.epsilon);
            gamma1 = @(tau) base.gamma1f(par, hpolicy(tau), taufunc_nu(tau), par.epsilon);
            delta = @(tau) base.deltaf(par, hpolicy(tau), taufunc_nu(tau), par.epsilon);
            labinc = @(tau) log(1-tau+par.X^(par.xi)*base.Finf(par,hpolicy(tau)));
        
            tax_choice = NaN(length(h_grid),1);
            nui = nu_grid(i);
            parfor j=1:length(h_grid)
                ObjectiveF = @(tau) -(par.omega*(h_grid(j)*log(gamma1_m(tau,h_grid(j)))+(1-h_grid(j))*log(gamma0_m(tau,h_grid(j))))+...
                    nui*(log(1-delta(tau))+par.beta*(log(delta(tau))+hpolicy(tau)*log(gamma1(tau))+(1-hpolicy(tau))*log(gamma0(tau)))+(1+par.beta)*labinc(tau))); %#ok<PFBNS>
                tax_choice(j) = fminbnd(ObjectiveF,-eps,1+eps,options);
            end
            tax_cell{i} = tax_choice;
        end
    end

    
    
end
end
    