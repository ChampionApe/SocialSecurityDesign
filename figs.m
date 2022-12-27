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
                    'Linewidth',2,'DisplayName','$\tau^u$');
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
                    'Linewidth',2,'DisplayName','$\delta_t^{u}(\tau_t^{c})$','Linestyle','--');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));                
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.srate(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$\delta_t^{u}$');
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
                    'Linewidth',2,'DisplayName','$h_t^{u}(\tau_t^{c})$','Linestyle','--');
         set(cu, 'MarkerEdgeColor',psettings.Mcolor(3,:),'MarkerFaceColor',psettings.Mcolor(3,:),...
                    'Color',psettings.color(3,:));
                
        u = plot(dates(1,psettings.start:psettings.end),sim_u.h(psettings.start:psettings.end),...
                    'Linewidth',2,'DisplayName','$h_t^{u}$');
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
            psettings.color = [[112/255, 128/255,144/255];[47/255,79/255,79/255];[0,0,128/255]];
            psettings.Mcolor=[[0,0,28/255];[0,0,28/255];[0,0,28/255]];
        end
    end
    
end
end