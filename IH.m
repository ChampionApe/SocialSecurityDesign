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
        % 2.1. Steady state policy function with EGM:
        function [sample,dif,maxdiff,x,tau_cell,tau_m_cell,h] = GradientFree_Policy(par,np,grids,hfunc,tau0)
        % 1: tau0 is a first guess on a policy function
        % tau(tau_{t-1},nu_t). If nan, guess constant function:
        if isnan(tau0)
            tau0 = @(tau,nu) 0.5-0.1.*tau;
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
        
        % Aux function - taufunction
        function taufunc = update_taufunc(tau_array,tau_m_array,tau,i)
            [tau_m,index] = sort(tau_m_array{i});
            func = griddedInterpolant(tau_m, tau_array{i}(index),'pchip','linear');
            taufunc = func(tau);
        end
        
        % Aux function - identify optimal policy along a grid given
        % continuation policy
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
                    ObjectiveF = @(tau) -(par.omega*(h_grid(j).*log(gamma1_m(tau,h_grid(j)))+(1-h_grid(j)).*log(gamma0_m(tau,h_grid(j))))+...
                        nui.*(log(1-delta(tau))+par.beta.*(log(delta(tau))+hpolicy(tau).*log(gamma1(tau))+(1-hpolicy(tau)).*log(gamma0(tau)))+(1+par.beta).*labinc(tau))); %#ok<PFBNS>
                    tax_choice(j) = fminbnd(ObjectiveF,-eps,1+eps,options);                
                end
                tax_cell{i} = tax_choice;
            end
        end
    
        % 2.2: Solve EGM-like given terminal policy function
        function [sol,par] = identify_policy_givenT(par, T_tau, T_taum, T_h)
            % 1: Solution structure
            sol = struct();
            [sol.tau, sol.tau_m, sol.h] = deal(cell(par.T,1),cell(par.T,1),cell(par.T,1));
            % 2: Add terminal period policy:
            sol.h{par.T} = num2cell(repmat(T_h,1,par.Ndata),1);
            sol.tau{par.T} = num2cell(repmat(T_tau,1,par.Ndata),1);
            sol.tau_m{par.T} = num2cell(repmat(T_taum,1,par.Ndata),1);
            % Create new grid of h to solve over:
            h = num2cell(repmat(SolveF.RefineTerminalGrid(par,par.h_upper,par.h_lower,2), 1, par.Ndata),1);
            for t = par.T-1:-1:1
                % 1. exogenous grid of h
                sol.h{t} = h; % for simplicity, keep grid-refinement from terminal period.
                % 2. Optimal tax in t, given policy function t+1 (use constant
                % epsilon)
                sol.tau{t} = SolveF.tau2cell(par,sol,t,par.epsilon,par.epsilon);
                % 3: Back out tau_m from grid of tau,h; vectorized, const.
                % epsilon and theta
                tau_m = base.BOT(par,vertcat(sol.h{t}{:}),vertcat(sol.tau{t}{:}),par.epsilon,par.theta);
                % Remove nan: (nb: vectorize this at some point)
                count = 0;
                for i=1:par.Ndata
                    sol.tau_m{t}{i} = tau_m(count+1:count+size(sol.h{t}{i},1));
                    count = count+size(sol.h{t}{i},1);
                    sol.tau{t}{i} = sol.tau{t}{i}(~isnan(sol.tau_m{t}{i}));
                    sol.h{t}{i} = sol.h{t}{i}(~isnan(sol.tau_m{t}{i}));
                    sol.tau_m{t}{i} = rmmissing(sol.tau_m{t}{i});
                end
            end
        end
        
        % 2.3: Solve EGM-like including terminal policy function. Combines
        % functions 2.1-2.2
        function [sol,par,maxdiff,x] = identifyPolicy(par, np, grids, hfunc, tau0)
            [~,~,maxdiff,x,tau_cell,tau_m_cell,h] = IH.GradientFree_Policy(par,np,grids,hfunc,tau0);
            [sol,par] = IH.identify_policy_givenT(par,tau_cell{1}, tau_m_cell{1}, h{1});
        end

        % 2.4: Static calibration with steady state assumption - EGM-like
        % solver:
        function [par,par_before,hfunc] = cal_ss_c(par, target, approx, np)
           par_before = par;
           par.alpha = 1-target.LSI;
           par.A = target.r/par.alpha;
           par.beta = base.cal_ss_beta(par,target);
           par.X = base.cal_ss_X(par,target);
           % Based on these parameters, approximate the labor equilibrium
           % function on a grid of tau_t, tau_{t+1}:
           hfunc = IH.h_EE_approximation(par,approx);
           np.n_nu = 1;
           grids = struct();
           grids.h = linspace(hfunc(1,0),hfunc(0,1),np.n_h);
           grids.nu = par.nu(target.nu_row,target.nu_col);
           grids.test = linspace(0.01, 0.99, np.n_test);
           % Solve for a policy
           taufunc = @(omega) target.pensiontax - IH.ssCalibration_aux(par, np, grids, hfunc, omega, target);
           %par.omega = fsolve(taufunc, 1.1205); % with target.pensiontax =
           %0.271 and h = 0.68
           %par.omega = fsolve(taufunc, 0.89685); % with target.pensiontax =
           %0.271 and h = 0.626
           par.omega = fsolve(taufunc, 0.804418);
        end
        
        % 2.5: Auxiliary function for 2.4
        function tau = ssCalibration_aux(par, np, grids, hfunc, omega,target)
            par.omega = omega;
            [~, ~, ~, ~, tau_cell, ~, h] = IH.GradientFree_Policy(par, np, grids, hfunc, nan);
            f = griddedInterpolant(h{1}, tau_cell{1}, 'pchip','linear');
            tau = f(target.h);
        end
        
        % 2.6: Steady state approximation simulation:
        
        
    % 3: Steady state approximation - VFI:
        % 3.1: Policy from grids
        function [sample, dif, maxdiff,x, tax, tax_m, h] = VFI_Policy(par,np,grids,tau0)
            % 1: tau0 is a first guess on a policy function
            % tau(tau_{t-1},nu_t). If nan, guess constant function:
            if isnan(tau0)
                tau0 = @(h,i) 0.5;
            end
            % 2: Data structure:
            maxdiff = NaN(np.max_iter,1); % difference between iterations
            [sample,dif] = deal(cell(np.max_iter,1)); % samples
            [tax,h] = deal(num2cell(zeros(length(grids.h),length(grids.nu)),1)');
            tax_m = deal(cell(np.n_nu,1));
            % 3: Iteration:
            for x=1:np.max_iter
                if x>1
                    tau0 = @(h,i) taufunc(h,i); % i refers to index in the nu_grid.
                end
                % Step 1: Given policy guess, find optimal tau:
                [tax,h] = IH.polobj_infinite_VFI(par,np,grids.nu,grids.h,tau0, tax, h,x);
                % Step 2: Check out feasible allocations:
                tau_m_vec = base.BOT(par,repmat(grids.h',np.n_nu,1),vertcat(tax{:}),par.epsilon,par.theta);
                % Step 3: Remove NaN (nb: no need for vectorization here)
                for i=1:length(tax)
                    tax_m{i} = tau_m_vec((i-1)*np.n_h+1:i*np.n_h);
                    tax_clean{i} = tax{i}(~isnan(tax_m{i}));
                    h_clean{i} = grids.h(~isnan(tax_m{i}))';
                    tax_m_clean{i} = rmmissing(tax_m{i});
                end
                
                % Step 4: Update the policy function: 
                taufunc = @(h,i) IH.update_taufunc_VFI(h_clean,tax_clean,h,i);
                
                % Step 5: Test convergence of function
                sample{x} = NaN(np.n_test, np.n_nu);
                dif{x} = NaN(np.n_test, np.n_nu);
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

        % 3.2: Political objective, infinite horizon, VFI:
        function [tax_cell, h_cell] = polobj_infinite_VFI(par,np,nu_grid,h_grid,taufunc, Ltax, Lh, x)
            [tax_cell, h_cell] = deal(cell(length(nu_grid),1));
            options = optimoptions('fmincon','ConstraintTolerance',np.tol_solver,'OptimalityTolerance',np.tol_solver);

            for i=1:length(nu_grid)
                taufunc_nu = @(h) taufunc(h,i);
                gamma0_m = @(tau,hi) base.gamma0f(par,hi,tau,par.epsilon);
                gamma1_m = @(tau,hi) base.gamma1f(par,hi,tau,par.epsilon);
                gamma0 = @(h) base.gamma0f(par, h, taufunc_nu(h), par.epsilon);
                gamma1 = @(h) base.gamma1f(par, h, taufunc_nu(h), par.epsilon);
                delta = @(h) base.deltaf(par, h, taufunc_nu(h), par.epsilon);
                labinc = @(tau,h) log(1-tau+par.X^(par.xi)*base.Finf(par,h));
       
                tax_choice = NaN(length(h_grid),1);
                h_choice = NaN(length(h_grid),1);
                nui = nu_grid(i);
                Ltaxi = Ltax{i};
                Lhi = Lh{i};
                parfor j=1:length(h_grid)
                    ObjectiveF = @(x) -(par.omega*(h_grid(j)*log(gamma1_m(x(1),h_grid(j)))+(1-h_grid(j))*log(gamma0_m(x(1),h_grid(j))))+...
                        nui*(log(1-delta(x(2)))+par.beta*(log(delta(x(2)))+x(2)*log(gamma1(x(2)))+(1-x(2))*log(gamma0(x(2))))+(1+par.beta)*labinc(x(1),x(2)))); %#ok<PFBNS>
                    lb = [eps, eps];
                    ub = [1-eps, 1-eps];
                    
                    nonlcon = @(x) IH.hConstraint(par, taufunc_nu, x(1), x(2));
                    sol = fmincon(ObjectiveF, [0.35,0.7], [], [], [],[], lb, ub, nonlcon, options);
%                    if x == 1
%                        sol = fmincon(ObjectiveF, [0.35,0.7], [], [], [],[], lb, ub, nonlcon, options);
%                    else
%                        sol = fmincon(ObjectiveF,[min(max(Ltaxi(j),eps),1-eps),min(max(Lhi(j),eps),1-eps)],[],[],[],[],lb,ub,nonlcon,options);
%                    end
                    tax_choice(j) = sol(1);
                    h_choice(j) = sol(2);
                end
                tax_cell{i} = tax_choice;
                h_cell{i} = h_choice;
            end
        
        end
    
        function [c,ceq] = hConstraint(par,taufunc_nu,tau,h)
            c = [];
            ceq = h-base.hfunc(par,h,tau,taufunc_nu(h),par.epsilon,par.theta);
        end
                
        function taufunc = update_taufunc_VFI(hgrid,tax,h,i)
            func = griddedInterpolant(hgrid{i}, tax{i},'pchip');
            taufunc = func(h);
        end
        
end
end
    