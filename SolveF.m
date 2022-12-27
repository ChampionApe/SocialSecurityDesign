classdef SolveF
methods(Static)
    %% 1: EGM finite horizon functions
    
    % 1.1: Create main struct with parameter values, exogenous population
    % data etc.
    function par = setup(setting)
        par = struct();
        par.A = 8;
        par.alpha = 0.47; % capital-income share
        par.omega = 0.8; % political weight of retirees
        par.xi = 0.35; % elasticity of labor supply
        par.X = 5; % measure of informal productivity. Requires X >1.
        par.beta = 0.326; % impatience parameter
        par.Gridpoints = 250; %number of gridpoints in 'solve'.
        par.Nsim = 100; % number of different h0 to simulate time-path.
        par.h0_sim = 80; % which index in the grid of h0 should we use to simulate.
        %par.TerminalAdhocAdjust=1;
        par.TerminalAdhocAdjust=(1+par.beta); % Set this to 1 to get the generic model. We might want to adjust upwards, to avoid a corner solution in terminal period.
        if setting == 'S'
            % Turn this on for some simple comparative statics of the model (for
            % some of the figures)rge:
            par.T=4;
            par.nu= 1.2*ones(par.T,1); % Turn on to illustrate figures.
        elseif setting == 'E'
            % Read data on population growth rates:
            [Rawdata, Names] = xlsread('Inputdata', '30year');
            PopulationIndex = strcmp(Names,'nu_t');
            YearIndex = strcmp(Names, 'Date');
            par.nu= rmmissing(Rawdata(:,PopulationIndex));
            par.T=length(par.nu); 
            par.years = rmmissing(Rawdata(:,YearIndex));
            par.Ndata = size(par.nu,2);
            par.dates = xlsread('Inputdata', 'PopulationData', 'A2:A37')';

            
            % Vectorized version:
            par.nu_v = par.nu';
            par.nu_v = par.nu_v(:);
            par.T_v = length(par.nu_v);
        end
    end
    
    % 1.2: Define upper/lower bounds:
    function [h_upper, h_lower] = h_bounds(par,tau_upper,tau_lower,tau_p_upper, tau_p_lower)
        h_upper = base.h_equi(par,tau_lower,tau_p_upper,par.epsilon,par.theta,0.5);
        h_lower = base.h_equi(par,tau_upper,tau_p_lower,par.epsilon,par.theta,0.5);
    end
    
    % 1.3: Add upper/lower bounds on variables to the parameter struct
    function par = bounds(par)
        [par.h_upper, par.h_lower] = SolveF.h_bounds(par,1-eps,eps,1-eps,eps);
    end
    
    % 1.4: Given parameters (par), identify policy functions:
    function [sol,par] = identify_policy(par)
       % 1: Solution structure
       sol = struct();
       [sol.tau, sol.tau_m, sol.h] = deal(cell(par.T,1),cell(par.T,1),cell(par.T,1));
       % 2: bounds:
       par = SolveF.bounds(par);
       % 3: Loop through time:
       for t=par.T
           % Define initial grids:
           [h,tau,tau_m] = SolveF.RefineTerminalGrid(par,par.h_upper,par.h_lower,2);
           
           % Create multidimensional array for vectorized operations:
           sol.h{t} = num2cell(repmat(h,1,par.Ndata),1);
           sol.tau{t} = num2cell(repmat(tau,1,par.Ndata),1);
           sol.tau_m{t} = num2cell(repmat(tau_m,1,par.Ndata),1);
       end
       for t = par.T-1:-1:1
           % 1. exogenous grid of h
           sol.h{t} = sol.h{par.T}; % for simplicity, keep grid-refinement from terminal period.
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

    % 1.5: Simulate path given policy functions and parameters            
    function sim = sim_u(par)
        sim = struct();
        sim.tau = base.tau_u(par,par.nu_v);
        sim.h = SolveF.h_universal(par,sim.tau);
        sim.srate = SolveF.compute_srate(par,sim,0);
    end
    function sim = sim_c(par,sol)
        sim = struct(); 
        h0 = linspace(par.h_lower,par.h_upper,par.Nsim)';
        sim.h = NaN(size(h0,1),par.T*par.Ndata);
        sim.tau = NaN(size(sim.h));
        % given h0, define PEE:
        for i=1:par.Ndata
            taxpolicy = griddedInterpolant(sol.h{1}{i},sol.tau{1}{i},'pchip','linear');
            sim.tau(:,i) = taxpolicy(h0);
        end
        for t=2:par.T
            for i=1:par.Ndata
                [tau_m,index] = sort(sol.tau_m{t}{i});
                taxpolicy_tau_m = griddedInterpolant(tau_m,sol.tau{t}{i}(index),'linear','none');
                h_tau_m = griddedInterpolant(tau_m,sol.h{t}{i}(index),'linear','none');
                sim.tau(:,(t-1)*par.Ndata+i) = taxpolicy_tau_m(sim.tau(:,(t-2)*par.Ndata+i));
                sim.h(:,(t-2)*par.Ndata+i) = h_tau_m(sim.tau(:,(t-2)*par.Ndata+i));
            end
        end
        % Terminal labor supply is defined as if universal+beveridgean:
        sim.h(:,(par.T-1)*par.Ndata+1:end) = base.h_equi(par,sim.tau(:,(par.T-1)*par.Ndata+1:end),1,0,0,0.5);
    end
    function sim = sim_c_unique(par,sol)
        sim_all = SolveF.sim_c(par,sol);
        sim = struct();
        sim.h = sim_all.h(par.h0_sim,:)';
        sim.tau = sim_all.tau(par.h0_sim,:)';
        sim.srate = SolveF.compute_srate(par,sim,par.epsilon);
    end
    function sim = sim_u_GivenTaxes(par,tau)
        sim = struct();
        sim.tau = tau;
        sim.h = SolveF.h_universal(par,tau);
        sim.srate = SolveF.compute_srate(par,sim,0);
    end
    function sim = sim_c_GivenTaxes(par,tau)
        sim = struct();
        sim.tau = tau;
        sim.h = NaN(size(tau));
        sim.h(1:end-1) = base.h_equi(par,tau(1:end-1),tau(2:end), par.epsilon,par.theta,zeros(size(tau(1:end-1)))+0.5);
        sim.h(end) = base.h_equi(par,tau(end),1,0,0,0.5); % terminal policy
        sim.srate = SolveF.compute_srate(par,sim,par.epsilon);
    end
    
    % 1.6: Approximate steady state policy with contributive pensions:
    function policy = ss_policy_c(par,nu)
        par.T = 6; % how many periods approximate ss?s
        par.Ndata = 1; 
        par.nu = nu*ones(par.T,1); % copy of nu data at steady state level
        policy = SolveF.identify_policy(par);
    end
    
    % Vectorized h for universal case:
    function h = h_universal(par,tau)
       h = NaN(size(tau));
       h(1:end-1) = base.h_equi(par, tau(1:end-1), tau(2:end), 0, par.theta, zeros(size(tau(1:end-1)))+0.5);
       h(end) = base.h_equi(par,tau(end),1,0,0,0.5);
    end
            
    % 1.7: Refine grid: 
    function [h,tau,tau_m] = RefineTerminalGrid(par, h_upper, h_lower,N)
        for j = 1:N
            if j~=1
                [h_upper,h_lower] = SolveF.h_bounds(par,1-eps, eps, max(tau(~isnan(tau_m))), min(tau(~isnan(tau_m))));
            end
            h = linspace(h_lower,h_upper,par.Gridpoints)';
            tau = base.tauT(par,h,par.epsilon);
            tau_m = base.BOT(par,h,tau,par.epsilon,par.theta);
            % remove nan:
            tau = tau(~isnan(tau_m));
            h = h(~isnan(tau_m));
            tau_m = rmmissing(tau_m);
        end
    end
    
    % 1.8: Collect optimal tax solution in cell:
    function tau = tau2cell(par,sol,t,epsilon,epsilon_p)
        tau = cell(1,par.Ndata);
        for i=1:par.Ndata
            tau{i} = base.tau_c(par, sol.h{t}{i}, sol.tau_m{t+1}{i}, sol.tau{t+1}{i}, sol.h{t+1}{i}, par.nu(t,i), epsilon, epsilon_p);
        end
        
    end
 
    % 1.9: Compute savings rate base on simulated time path
    function srate = compute_srate(par,sim,epsilon)
        gamma0 = base.gamma0f(par,sim.h(1:(end-1)),sim.tau(2:end),epsilon);
        gamma1 = base.gamma1f(par,sim.h(1:(end-1)),sim.tau(2:end),epsilon);
        srate = base.deltaf(par,sim.h(1:(end-1)),gamma0,gamma1);
    end
    
    
    %% 2: Calibration with EGM finite horizon, contributive
    % 2.1: Default calibration targets:
    function t = target()
        t = struct();
        t.LSI = 0.5;
        t.r = 4.3; % 30-year cumulated interest rate.
        t.srate = 0.206;
        t.pensiontax = 0.271;
        t.h = 0.5777; % 1-informality rate
        t.nu = 1.2120; % 1.2120 corresponding to level in 2010 when data starts in 1950.
        t.nu_row = 3; % which row in par.nu is the target in.
        t.nu_col = 1; % which column in par.nu is the target in.
        t.xi = 0.35;
    end

    % 2.2: Approximate steady state and calibrate with contributive
    % pensions:
    function [par, par_before] = cal_ss_c(par,target)
        par_before = par;
        par.alpha = 1-target.LSI;
        par.A = target.r/par.alpha;
        % Back out beta from savings rate:
        par.beta = base.cal_ss_beta(par,target);
        % Solve for X
        par.X = base.cal_ss_X(par,target);
        % Solve for omega by targeting the simumlated pension rate:
        taufunc = @(omega) target.pensiontax-base.cal_ss_omega(par,target,omega);
        par.omega = fsolve(taufunc, par.omega);
    end
    
    % 2.3: Dynamic calibration, finite horizon EGM, contributive:
    function [par, par_before] = cal_dyn_c(par,target)
        par_before = par;
        par.alpha = 1-target.LSI;
        par.A = target.r/par.alpha;
        % Simulation calibration:
        par_temp = par;
        [par_temp.nu,par_temp.Ndata] = deal(par.nu(:,target.nu_col), 1);
        
        % Function evaluating to zeros when targets are met:
        vec_func = @(x) SolveF.dyn_aux(par_temp,target,x);
        [xsol,~,exitflag] = fsolve(vec_func, [par.beta, par.X, par.omega]);
        assert(exitflag>0, 'Calibration failed');
        % Update parameters:
        [par.beta, par.X, par.omega] = deal(xsol(1), xsol(2), xsol(3));
    end
    
    % Auxiliary function used for dynamic calibration:
    function zeros = dyn_aux(par,target,x)
        zeros = NaN(3,1); % solve for three parameters
        [par.beta, par.X, par.omega] = deal(x(1),x(2),x(3));
        %simulate model:
        [sol_temp,par] = SolveF.identify_policy(par);
        sim_temp = SolveF.sim_c_unique(par,sol_temp);
        % relevant variables:
        [h1980, tax2010] = deal(sim_temp.h(target.nu_row-1), sim_temp.tau(target.nu_row));
        gamma0 = base.gamma0f(par,h1980,tax2010,par.epsilon);
        gamma1 = base.gamma1f(par,h1980,tax2010,par.epsilon);
        delta2010 = base.deltaf(par,h1980,gamma0,gamma1);
        zeros(1) = h1980-target.h;
        zeros(2) = tax2010-target.pensiontax;
        zeros(3) = delta2010-target.srate;
    end

    
    
    
end
end
