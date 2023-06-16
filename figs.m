classdef figs
methods(Static)
    
    % Plot the policy function:
    function [] = policyfunction(dates,sol,t,print,prefix)
        
        fig=figure('Name', sprintf('%s_policyfunction', prefix));
        hold('on')
        
        % rearrange 'dates' as the solution struct is:
        dates = reshape(dates, [length(sol.tau), length(sol.tau{1})])';
        % Identify the correct entry for year 't':
        [a,b] = find(dates==t);
        
        hgrid = sol.h{a}{b};
        taugrid = sol.tau{a}{b};
        
        tau = plot(hgrid,taugrid,...
                       'Linewidth', 2, 'DisplayName', '$\tau_t$');
                         
                      set(tau, 'MarkerEdgeColor', [0 0 28/255], 'MarkerFaceColor', [0 0 28/255],...
                        'Color', [112/255 128/255 144/255]); 
                   
         xlabel('$h_{t-1}$')
         ylabel('$\tau_t$')
         legend('Location','northwest');
         box('on');
         grid on;
         
         if print==1
              base.printfig(fig);
         end
         
    end
    
    function [] = plot_tax(par,dates,sim,print,prefix)
        
        fig=figure('Name', sprintf('%s_tax', prefix));
        hold('on')
        
        par.nu = par.nu_v;
        par.T = par.T_v;

        % Start/end plots at indices:
        startplot = 7;
        endplot=par.T_v-9;
        dates=dates(1,startplot:endplot);

        taxpath = sim.tau(startplot:endplot);
        
        tau = plot(dates,taxpath,...
                       'Linewidth', 2, 'DisplayName', '$\tau_t^c$');
                         
                      set(tau, 'MarkerEdgeColor', [0 0 28/255], 'MarkerFaceColor', [0 0 28/255],...
                        'Color', [112/255 128/255 144/255]); 
                   
         xlabel('Years')
         ylabel('$\tau_t$')
         legend('Location','northwest');
         box('on');
         grid on;
         
         if print==1
              base.printfig(fig);
         end
         
    end
    
    function [] = w2r_ratio(dates,par,psettings,print,prefix)
        
        fig=figure('Name', sprintf('%s_w2r', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        v = plot(dates(1,psettings.start:psettings.end),par.nu_v(psettings.start:psettings.end),...
                    'Linewidth',2,'DefaultAxesFontSize', psettings.FontSize);
        set(v, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:)); 
        set(gca,'FontSize',psettings.FontSize)
         xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\nu_t$','FontSize', psettings.LabelFontSize)
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end

    function [] = comparetaxes(dates,psettings,sim_c,sim_u,print,prefix)
        
        fig=figure('Name', sprintf('%s_tax', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau^c$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));
        
        u = plot(dates(1,psettings.start:psettings.end),sim_u.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau^u$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\tau_t$','FontSize', psettings.LabelFontSize)
        legend('Location','northwest','FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end
    
    function [] = comparesavings(dates,psettings,sim_c,sim_u,sim_counterfactual,print,prefix)
        
        fig=figure('Name', sprintf('%s_srate', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t^c$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));

        cu = plot(dates(1,psettings.start:psettings.end),sim_counterfactual.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t^{u}(\tau_t^{c})$','Linestyle','-.');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));                
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t^{u}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
                
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\delta_t$','FontSize', psettings.LabelFontSize)
        legend('Position',[0.65 0.7 0.18 0.18],'FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end

    function [] = comparelabor(dates,psettings,sim_c,sim_u,sim_counterfactual,print,prefix)
        
        fig=figure('Name', sprintf('%s_h', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t^{c}$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));
                
        cu = plot(dates(1,psettings.start:psettings.end),sim_counterfactual.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t^{u}(\tau_t^{c})$','Linestyle','-.');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t^{u}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
                
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$h_t$','FontSize', psettings.LabelFontSize)
        legend('Position',[0.65 0.7 0.18 0.18],'FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end
    
    function [] = SurfPlot(psettings,epsgrid,thetagrid,variable,name,print,prefix)
       fig = figure('Name',sprintf('%s_grids_%s', prefix,name));
       hold('on')
       psettings = figs.def_psettings(psettings);
       [X,Y] = meshgrid(thetagrid,epsgrid);
       surf(X,Y,variable,...
            'LineStyle', ':', 'FaceAlpha', 0.75);
       colormap bone;
       xlabel('$\theta$','FontSize', psettings.LabelFontSize);
       ylabel('$\epsilon$','FontSize', psettings.LabelFontSize);
       if strcmp(name, 'tax')
           view([-111.918 34.477]);
           zlabel('$\tau$', 'FontSize', psettings.LabelFontSize);
       elseif strcmp(name,'srate')
           view([-56.3439 20.4059]);
           zlabel('$\delta$', 'Fontsize', psettings.LabelFontSize);
       elseif strcmp(name,'h')
           view([-109.0435 37.6750]);
           zlabel('$h$', 'Fontsize', psettings.LabelFontSize);
       end
       box('on');
       grid on;
       if print==1
           base.printfig(fig);
       end
    end

    function [] = comparetaxes_US(dates,psettings,sim_c,sim_u,print,prefix)
        
        fig=figure('Name', sprintf('%s_tax', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:)); 
        u = plot(dates(1,psettings.start:psettings.end),sim_u.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau^{covid}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\tau_t$','FontSize', psettings.LabelFontSize)
        legend('Location','northwest','FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end
    
    function [] = comparesavings_US(dates,psettings,sim_c,sim_u,sim_counterfactual,print,prefix)
        
        fig=figure('Name', sprintf('%s_srate', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));

        cu = plot(dates(1,psettings.start:psettings.end),sim_counterfactual.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t(\tau_t^{covid})$','Linestyle','-.');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));                
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t^{covid}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
                
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\delta_t$','FontSize', psettings.LabelFontSize)
        legend('Position',[0.65 0.7 0.18 0.18],'FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end

    function [] = comparelabor_US(dates,psettings,sim_c,sim_u,sim_counterfactual,print,prefix)
        
        fig=figure('Name', sprintf('%s_h', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));
                
        cu = plot(dates(1,psettings.start:psettings.end),sim_counterfactual.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t(\tau_t^{covid})$','Linestyle','-.');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t^{covid}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
                
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$h_t$','FontSize', psettings.LabelFontSize)
        legend('Position',[0.65 0.7 0.18 0.18],'FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end

    function [] = policyConvergence(dates, taxFunctions,print,prefix)
        
    end
    
    function [] = comparetaxes_Convergence(dates,psettings,sims,grids,print,prefix)
        
        fig=figure('Name', sprintf('%s_taxConvergence', prefix));
        hold('on')
        set(gca, 'ColorOrder', psettings.color);
        for i=1:length(grids.eps)
            plot(dates(1,psettings.start:psettings.end), sims(psettings.start:psettings.end,i),...
                'Linewidth', 2, 'LineStyle', psettings.lineStyles(i), 'DisplayName', sprintf('$\\epsilon=$ %0.2g', grids.eps(i)));
        end
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\tau_t$','FontSize', psettings.LabelFontSize)
        legend('Location','northwest','FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
        if print==1
              base.printfig(fig);
         end
    end
   
    function [] = comparetaxes_IH_FH(dates,psettings,sim_c,sim_u,print,prefix)
        
        fig=figure('Name', sprintf('%s_tax', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);
        
        c = plot(dates(1,psettings.start:psettings.end),sim_c.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau^{FH}$');
        set(c, 'MarkerEdgeColor', psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
                    'Color', psettings.color(1,:));
        
        u = plot(dates(1,psettings.start:psettings.end),sim_u.tau(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\tau^{IH}$','Linestyle','--');
         set(u, 'MarkerEdgeColor',psettings.Mcolor(2,:),'MarkerFaceColor',psettings.Mcolor(2,:),...
                    'Color',psettings.color(2,:));
        set(gca,'FontSize',psettings.FontSize);
        xlabel('Years','FontSize', psettings.LabelFontSize)
        ylabel('$\tau_t$','FontSize', psettings.LabelFontSize)
        legend('Location','northwest','FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end
    end
    
    function [] = SSpolicies(psettings, solTau, solTau_m, threshold, grids, print, prefix)
        
        fig=figure('Name', sprintf('%s_SS_policy', prefix));
        hold('on')
        % Add default settings to psettings:
        psettings = figs.def_psettings(psettings);

        z = plot(solTau_m{threshold-1}, solTau{threshold-1},...
            'Linewidth', 2, 'Linestyle','-.','DisplayName', sprintf('$\\nu=$ %.2f', grids.nu(threshold-1)));
        set(z, 'MarkerEdgeColor',psettings.Mcolor(3,:), 'MarkerFaceColor', psettings.Mcolor(3,:),...
            'Color', psettings.color(3,:));
        
        x = plot(solTau_m{threshold}, solTau{threshold},...
            'Linewidth', 2, 'DisplayName', sprintf('$\\nu=$ %.2f', grids.nu(threshold)));
        set(x, 'MarkerEdgeColor',psettings.Mcolor(1,:), 'MarkerFaceColor', psettings.Mcolor(1,:),...
            'Color', psettings.color(1,:));

        y = plot(solTau_m{threshold+1}, solTau{threshold+1},...
            'Linewidth', 2, 'Linestyle','--','DisplayName', sprintf('$\\nu=$ %.2f', grids.nu(threshold+1)));
        set(y, 'MarkerEdgeColor',psettings.Mcolor(2,:), 'MarkerFaceColor', psettings.Mcolor(2,:),...
            'Color', psettings.color(2,:));
        
        set(gca,'FontSize',psettings.FontSize);
        xlabel('$\tau_{t-1}$','FontSize', psettings.LabelFontSize)
        ylabel('$\tau_t$','FontSize', psettings.LabelFontSize)
        legend('Location','southwest','FontSize',psettings.LabelFontSize);
        box('on');
        grid on;
         if print==1
              base.printfig(fig);
         end

    end
    
    function psettings = def_psettings(psettings)
        if ~isfield(psettings,'start')
            psettings.start = 7;
        end
        if ~isfield(psettings,'end')
            psettings.end= 27;
        end
        if ~isfield(psettings,'FontSize')
            psettings.FontSize= 15;
        end
        if ~isfield(psettings,'LabelFontSize')
            psettings.LabelFontSize= 20;
        end
        if ~isfield(psettings,'color')
            psettings.color = [[0.12156862745098039, 0.4666666666666667, 0.7058823529411765];...
                                [1.0, 0.4980392156862745, 0.054901960784313725];...
                                [0.17254901960784313, 0.6274509803921569, 0.17254901960784313];...
                                [0.8392156862745098, 0.15294117647058825, 0.1568627450980392];...
                                [0.5803921568627451, 0.403921568627451, 0.7411764705882353]];
            psettings.Mcolor=[[0,0,28/255];[0,0,28/255];[0,0,28/255]];
            psettings.lineStyles = ["-"; "--"; ":"; "-.", ;"-"];
        end
    end
    
end
end