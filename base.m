classdef base
methods(Static)

    % 1: EE FUNCTIONS 
    
    % Informal labor function:
    function F = Finf(par,h)
        F = par.X*(par.xi./(1+par.xi)).*(1-h.^(1+1/par.xi));
    end
    
    % Parameter functions:
    % g(h_t) (balancing budget condition)
    function gh = ghf(h,epsilon_p)
       gh = 1./(1+epsilon_p.*(2.*h-1));
    end
    
    % Pension multipliers:
    function gamma0 = gamma0f(par,h,tau_p,epsilon_p)
       gamma0 = 1+((1-par.alpha)./par.alpha).*tau_p.*base.ghf(h,epsilon_p).*(1-epsilon_p); 
    end
    function gamma1 = gamma1f(par,h,tau_p,epsilon_p)
        gamma1 = 1+((1-par.alpha)./par.alpha).*tau_p.*base.ghf(h,epsilon_p).*(1+epsilon_p);
    end
    
    % Savings rate:
    function delta = deltaf(par,h,gamma0,gamma1)
        delta = par.beta.*(gamma1.*(1-h)+gamma0.*h)./(par.beta.*(gamma1.*(1-h)+gamma0.*h)+gamma1.*gamma0);
    end
    
    % labor supply, auxiliary functions:
        % Common component for all systems:
        function param = hp_1(par,tau)
            param = (1-tau)./(par.X.^(1+par.xi));
        end
        % all systems except universal use this one:
        function param = hp_2(par,h,tau)
           param = par.beta.*h.*(base.hp_1(par,tau)+base.Finf(par,h)./par.X); 
        end
        
        % Gamma function:
        function param = hp_3(par,h,tau_p,epsilon_p,theta_p)
           gamma0 = base.gamma0f(par,h,tau_p,epsilon_p);
           gamma1 = base.gamma1f(par,h,tau_p,epsilon_p);
           param = (1-base.deltaf(par,h,gamma0,gamma1)).*(log(gamma1./gamma0)+theta_p.*((gamma1-1)./gamma1+((1-h)./h).*(gamma0-1)./gamma0));
        end

    % Labor equilibrium functions: 
    function chih = hfunc(par,h,tau,tau_p,epsilon_p,theta_p)
        if (theta_p == 0) && (epsilon_p == 0)
            chih = base.hp_1(par,tau).^(par.xi/(1+par.xi));
        else
            chih = (base.hp_1(par,tau)+base.hp_2(par,h,tau).*base.hp_3(par,h,tau_p,epsilon_p,theta_p)).^(par.xi/(1+par.xi));
        end
    end
    
    function h = h_equi(par,tau,tau_p,epsilon_p,theta_p,h0)
        if (theta_p == 0) && (epsilon_p == 0)
            h = base.hp_1(par,tau).^(par.xi/(1+par.xi));
        else
            h = fsolve(@(x) base.hfunc(par,x,tau,tau_p,epsilon_p,theta_p)-x,h0);
        end
    end
    
    % 2: TAX FUNCTIONS
    
    % 2.1: Terminal Tax functions
    function tau = tauT(par,h_m,epsilon)
        if epsilon == 0
            tau = (h_m./(par.TerminalAdhocAdjust.*par.nu(par.T)+par.omega.*h_m)).*(par.omega.*(1+par.xi.*par.X.^(1+par.xi))-par.TerminalAdhocAdjust.*par.nu(par.T).*(par.alpha./(1-par.alpha)));
        else
            gamma0 = @(x) base.gamma0f(par,h_m,x,epsilon); 
            gamma1 = @(x) base.gamma1f(par,h_m,x,epsilon);
            heq = @(x) x./(1-x+par.xi*par.X^(1+par.xi))-(par.omega./(par.nu(par.T).*(1+par.beta)))*(h_m.*(gamma1(x)-1)./gamma1(x)+(1-h_m).*(gamma0(x)-1)./gamma0(x));
            if isscalar(h_m)
                tau = fsolve(heq, 0.5);
            else
                tau = fsolve(heq, zeros(size(h_m))+0.5);
            end
        end
        tau = min(1,max(0,tau));
    end
        
    % 2.2: Tax function, contributive: h_grid, solution_p_Ltax,
    % solution_p_tax, solution_p_h, should all be vectors. nu,
    % epsilon,epsilon_p are scalars. 
    % Returns a vector of solutions corresponding to the h_grid.
    function tau = tau_c(par,h_grid,solution_p_Ltax, solution_p_tax, solution_p_h,nu,epsilon,epsilon_p)
        % 1: Interpolants:  This first section should only be included if
        % epsilon/theta can change over time.
        % 1.1: If the system in t+1 is universal, the policy function is
        % independent of current tax choices. In this case:
%        if isscalar(solution_p_tax)
%            taxpolicy = @(tau) solution_p_tax;
%            hpolicy = @(tau) solution_p_h;
%        elseif (sum(solution_p_tax == solution_p_tax(1)) == length(solution_p_tax))
%            taxpolicy = @(tau) solution_p_tax(1);
%            hpolicy = @(tau) solution_p_h(1);
%        % 1.2: If the system in t+1 is contributive, create a gridded
%        % interpolant for the policy and labor functions:
%        else
        % 1: Sort period t+1's lagged tax rate, tax rate, and labor supply:
        [Ltax, Index] = sort(solution_p_Ltax);
        tax = solution_p_tax(Index);
        h = solution_p_h(Index);
        % 2: Interpolation PEE tax and labor function:
        taxpolicy = griddedInterpolant(Ltax,tax,'pchip','none');
        hpolicy = griddedInterpolant(Ltax,h,'pchip','none');
%        end
        
        % 2: Auxiliary functions for the political objective:
        gamma0_minusone = @(tau,hi) base.gamma0f(par, hi, tau, epsilon);
        gamma1_minusone = @(tau,hi) base.gamma1f(par, hi, tau, epsilon);
        gamma0 = @(tau) base.gamma0f(par, hpolicy(tau), taxpolicy(tau), epsilon_p);
        gamma1 = @(tau) base.gamma1f(par, hpolicy(tau), taxpolicy(tau), epsilon_p);
        delta = @(tau) base.deltaf(par, hpolicy(tau), gamma0(tau), gamma1(tau));
        labinc = @(tau) log(1-tau+par.X^(par.xi)*base.Finf(par,hpolicy(tau)));
        
        % 3: Objective/solution:
        tau = NaN(length(h_grid),1);
        
        parfor j=1:length(h_grid)
            ObjectiveF = @(tau) -(par.omega.*(h_grid(j).*log(gamma1_minusone(tau,h_grid(j)))+(1-h_grid(j)).*log(gamma0_minusone(tau,h_grid(j))))+...
                nu.*(log(1-delta(tau))+par.beta.*(log(delta(tau))+hpolicy(tau).*log(gamma1(tau))+(1-hpolicy(tau)).*log(gamma0(tau)))+(1+par.beta).*labinc(tau))); %#ok<PFBNS>            
            tau(j) = fminbnd(ObjectiveF,-eps,1+eps);
        end
        
    end
        
    % 2.3: Tax function, universal; nu can be a vector of a scalar
    function tau = tau_u(par,nu)
       tau = min(1, max(0, (1./(par.omega+nu.*(1+par.beta))).*(par.omega.*(1+par.xi.*par.X.^(1+par.xi))-(par.alpha/(1-par.alpha)).*nu.*(1+par.beta)))); 
    end
    
    % 2.4: Tax function, inferring the solution based on values on epsilon:
    function tau = tauf(par,h_grid,solution_p_Ltax, solution_p_tax, solution_p_h,nu,epsilon,epsilon_p)
        if epsilon == 0
            tau = tau_u(par,nu);
        else
            tau = tau_c(par,h_grid,solution_p_Ltax, solution_p_tax, solution_p_h,nu,epsilon,epsilon_p);
        end
    end
    
    % 2.5: Back out lagged tax rate from tau and h using labor eq.:
    function tau = BOT(par,h,tau_p,epsilon_p,theta_p)
        Gamma = base.hp_3(par,h,tau_p,epsilon_p,theta_p);
        tau = 1-((par.X.^(1+par.xi))./(1+par.beta.*h.*Gamma)).*(h.^((1+par.xi)./par.xi)-par.beta.*h.*Gamma.*(par.xi./(1+par.xi)).*(1-h.^(1+1/par.xi)));
        tau(tau>1 | tau<0)=NaN; % Only return value if in [0,1].
    end

    % 3: Calibration:
    % 3.1: Steady state calibration functions:
    function beta = cal_ss_beta(par,target)
        gamma0 = base.gamma0f(par,target.h,target.pensiontax,par.epsilon);
        gamma1 = base.gamma1f(par,target.h,target.pensiontax,par.epsilon);
        beta = (target.srate/(1-target.srate))*(gamma0*gamma1)/(gamma1*(1-target.h)+gamma0*target.h);
    end
    function X = cal_ss_X(par,target)
        Gamma = base.hp_3(par, target.h,target.pensiontax,par.epsilon,par.theta);
        X = ((1-target.pensiontax)*(1+par.beta*target.h*Gamma)/(target.h^((1+par.xi)/par.xi)-par.beta*target.h*(par.xi/(1+par.xi))*(1-target.h^(1+1/par.xi))*Gamma))^(1/(1+par.xi));
    end
    function tau = cal_ss_omega(par,target,omega)
        par.omega = omega;
        soltemp = SolveF.ss_policy_c(par,target.nu);
        tax_policy = griddedInterpolant(soltemp.h{1}{1},soltemp.tau{1}{1},'pchip','linear');
        tau = tax_policy(target.h);
    end
    
    % Figure, layout, basic functions:
    function x = nonlinspace(lo,hi,n,phi)
        % recursively constructs an unequally spaced grid.
        % phi > 1 -> more mass at the lower end of the grid.
        % lo can be a vector (x then becomes a matrix).
        x      = NaN(n,length(lo));
        x(1,:) = lo;
        for i = 2:n
            x(i,:) = x(i-1,:) + (hi'-x(i-1,:))./((n-i+1)^phi);
        end
    end

    
    function [] = layout()    
        % set layout parameters
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaultTextInterpreter','latex');
        set(groot, 'defaultAxesFontSize', 12); 
        
    end
    
    function [] = printfig(figin)
        fig = figure(figin);
        fig.PaperUnits = 'centimeters';   
        fig.PaperPositionMode = 'manual';
        fig.PaperPosition = [0 0 16 12];
        fig.PaperSize = [16 12];

        filename = ['figs\' get(fig,'name') ''];            
        print('-dpdf',['' filename '.pdf']);

    end
   
end
end
