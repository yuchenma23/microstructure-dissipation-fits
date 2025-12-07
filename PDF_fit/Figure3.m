%%  LOG–NORMAL-MIXTURE FITS — three data sets, three panels each
% -------------------------------------------------------------------------
clear; close all; clc

load ../Dissipation_data_collect/eps_bbtre.mat          
load ../Dissipation_data_collect/eps_natre.mat

fontsize = 20;                                    % master font size


% ---------- 1) fit each data set -----------------------------------------
D.aby  = Fit_Mixture( epsilons_bbtre_0_1000);   
D.int  = Fit_Mixture( epsilons_bbtre_1000_2000 ); % lower interior (1000–2000 m)
D.nat  = Fit_Mixture( epsilons_natre_500_2000 );  % NATRE all depths

% ---------- 2) build the figure (now 3 × 3) ------------------------------
figure(1); clf
set(gcf,'color','w','Position',[1471 -111 1294 725])   % taller window

dataOrder = {'aby','int','nat'};                    % plotting order
rowTitles = {'Near-Bottom Layer (BBTRE)', ...
             'Lower Interior Layer (BBTRE)', ...
             'Upper Interior Layer (NATRE)'};
panelTag  = 'abcdefghi';                            % a) … i)

for r = 1:3
    S = D.(dataOrder{r});
    base = (r-1)*3;   % panel index offset

    % ---- (1) CDF --------------------------------------------------------
    subplot(3,3,base+1)
    plot(S.x , S.ecdf ,'k','LineWidth',3); hold on
    plot(S.xfit , S.cdf_mix ,'--','Color','red','LineWidth',3)
    text(0.02,0.95,sprintf('%c)',panelTag(base+1)),'Units','normalized',...
         'FontSize',fontsize)
    if r==3
        xlabel('$\log_{10}(\varepsilon \, [\mathrm{W}/\mathrm{kg}])) $','Interpreter','latex');
    end
    xlim([-13 -6]); ylim([0 1])
    if r==1
        legend('DATA','Mixture fit','Location','best','Interpreter','latex'); 
    end
    box off
    set(gca,'FontSize',fontsize,'FontName','times')
    ylabel('CDF')

    % ---- (2) PDF (linear-y) --------------------------------------------
    subplot(3,3,base+2)
       z = log10(S.raw);
    edges = linspace(-13, -6, 61);  % Nbin=60
    histogram(z, 'Normalization','pdf','BinEdges', edges, ...
          'EdgeColor',[.7 .7 .7], 'FaceColor',[.85 .85 .85]); hold on
    plot(S.xfit(2:end) , S.pdf_mix, '--', 'Color', 'red', 'LineWidth', 3)
    plot(S.xfit(2:end) , S.pdf_comp(1,:), ':', 'Color', [0 0.45 0.74], 'LineWidth', 2)
    plot(S.xfit(2:end) , S.pdf_comp(2,:), ':', 'Color', [0.00 0.50 0.00], 'LineWidth', 2)
    xline(log10(mean(S.raw)) ,'k','LineWidth',2)
    xline(log10(S.mean_mix) ,'Color','red','LineWidth',2)
    text(0.02,0.95,sprintf('%c)',panelTag(base+2)),'Units','normalized',...
         'FontSize',fontsize)
    if r==3
        xlabel('$\log_{10}(\varepsilon \, [\mathrm{W}/\mathrm{kg}])) $','Interpreter','latex');
    end
    if r==2
        legend({'Histogram','Mixture fit','Comp 1','Comp 2', ...
            '$\overline{\varepsilon}$','$\overline{\varepsilon}_{\mathrm{fit}}$'}, ...
            'Location','best','Interpreter','latex','FontSize',15)
    end
    ylabel('PDF')
    xlim([-13 -6]); ylim([1e-4 1])
    box off
    set(gca,'FontSize',fontsize,'FontName','times')

    % ---- (3) PDF (log-y) -----------------------------------------------
    subplot(3,3,base+3)
    z = log10(S.raw);
    edges = linspace(-13, -6, 61);  % Nbin=60
    histogram(z, 'Normalization','pdf','BinEdges', edges, ...
          'EdgeColor',[.7 .7 .7], 'FaceColor',[.85 .85 .85]);
    
    hold on

    plot(S.xfit(2:end) , S.pdf_mix, '--', 'Color', 'red', 'LineWidth', 3)
    plot(S.xfit(2:end) , S.pdf_comp(1,:), ':', 'Color', [0 0.45 0.74], 'LineWidth', 2)
    plot(S.xfit(2:end) , S.pdf_comp(2,:), ':', 'Color', [0.00 0.50 0.00], 'LineWidth', 2)
 %   xline(log10(mean(S.raw)) ,'k','LineWidth',3)
 %   xline(log10(S.mean_mix) ,'Color',[0.85 0.33 0.10],'LineWidth',3)
    text(0.02,0.95,sprintf('%c)',panelTag(base+3)),'Units','normalized',...
         'FontSize',fontsize)
    if r==3
        xlabel('$\log_{10}(\varepsilon \, [\mathrm{W}/\mathrm{kg}])) $','Interpreter','latex');
    end
    ylabel('PDF')
        yticks = -6:2:0;
    ytick_vals = 10.^yticks;
    set(gca, 'YTick', ytick_vals, ...
             'YTickLabel', arrayfun(@(x) sprintf('$10^{%d}$', x), yticks, 'UniformOutput', false), ...
             'TickLabelInterpreter','latex');
    xlim([-13 -6]); ylim([1e-4 5]); set(gca,'YScale','log')
    box off
    set(gca,'FontSize',fontsize,'FontName','times')
end

% ---------- 3) row-level titles ------------------------------------------
yPos = [0.93 0.63 0.34];                       % vertical placements
for k = 1:3
   annotation('textbox',[0.13 yPos(k) 0.74 0.04],...
      'String',rowTitles{k},...
      'EdgeColor','none','HorizontalAlignment','center',...
      'FontSize',fontsize+2,'FontWeight','bold','FontName','latex');
end

% ---------- 5) Save output -----------------------------------------------
print(gcf,'-djpeg','-r200','Figure3.jpg')


