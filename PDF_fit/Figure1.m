%%  LOG–SKEW-NORMAL FITS with MLE, three data sets, one plot
% -------------------------------------------------------------------------
clear;  close all;  clc
load ../Dissipation_data_collect/eps_bbtre.mat          
load ../Dissipation_data_collect/eps_natre.mat                  

fontsize = 20;                                    % master font


% ---------- 1) fit each depth range --------------------------------------
D.aby  = Fit_LSN( epsilons_bbtre_0_1000   );  % abyss     0–1000 m
D.int  = Fit_LSN( epsilons_bbtre_1000_2000  ); % interior 1000–2000 m
D.nat  = Fit_LSN( epsilons_natre_500_2000 );                 % NATRE all depths

% ---------- 2) build the figure (identical style, 3×3) -------------------
figure(1); clf
set(gcf,'color','w','Position',[1471 -111 1294 725])   % taller window

dataOrder = {'aby','int','nat'};                    % plotting sequence
rowTitles = {'Near-Bottom Layer (BBTRE)','Lower Interior Layer (BBTRE)','Upper Interior Layer (NATRE)'};
panelTag  = 'abcdefghi';                            % a)--> i)

for r = 1:3
    S = D.(dataOrder{r});
    base = (r-1)*3;          % column offset for this row (CDF | PDF | PDF-log)

    % --- (1) CDF ----------------------------------------------------------
    subplot(3,3,base+1)
    plot(S.x , S.ecdf ,'k','LineWidth',3); hold on
    plot(S.xfit , S.cdf_fit,'r--','LineWidth',3)
    text(0.02,0.95,sprintf('%c)',panelTag(base+1)),'Units','normalized',...
         'FontSize',fontsize)
    xlim([-13 -6]); ylim([0 1])
    if r==2, legend('DATA','LSN fit','Location','best','Interpreter','latex'); end
    if r==3, xlabel('$\log_{10}(\varepsilon\, [\mathrm{W}/\mathrm{kg}])$',...
                    'Interpreter','latex'); end
    ylabel('CDF')

    box off;
    set(gca,'Fontsize',fontsize,'fontname','times')

    % --- (2) PDF (linear-y) ----------------------------------------------
    subplot(3,3,base+2)
    z = log10(S.raw);
    edges = linspace(-13, -6, 61);  % Nbin=60
    histogram(z, 'Normalization','pdf','BinEdges', edges, ...
          'EdgeColor',[.7 .7 .7], 'FaceColor',[.85 .85 .85]);
    hold on
    plot(S.xfit(2:end) , S.pdf_fit,'r--','LineWidth',3)
    xline(log10(mean(S.raw)),'k','LineWidth',2)
    xline(log10(S.mean_fit) ,'r','LineWidth',2)
    text(0.02,0.95,sprintf('%c)',panelTag(base+2)),'Units','normalized',...
         'FontSize',fontsize)
    xlim([-13 -6]); ylim([1e-4 1])
    if r==2, legend('DATA','LSN fit','$\overline{\varepsilon}$',...
                    '$\overline{\varepsilon}_{\mathrm{fit}}$','Location','best',...
                    'Interpreter','latex'); end
    if r==3, xlabel('$\log_{10}(\varepsilon\, [\mathrm{W}/\mathrm{kg}])$',...
                    'Interpreter','latex'); end
    ylabel('PDF')
    box off;
    set(gca,'Fontsize',fontsize,'fontname','times')

    % --- (3) PDF (log-y) --------------------------------------------------
    subplot(3,3,base+3)
    z = log10(S.raw);
    edges = linspace(-13, -6, 61);  % Nbin=60
    histogram(z, 'Normalization','pdf','BinEdges', edges, ...
          'EdgeColor',[.7 .7 .7], 'FaceColor',[.85 .85 .85]);
    hold on
    plot(S.xfit(2:end) , S.pdf_fit,'r--','LineWidth',3)
    text(0.02,0.95,sprintf('%c)',panelTag(base+3)),'Units','normalized',...
         'FontSize',fontsize)
    xlim([-13 -6]); ylim([1e-4 5]); set(gca,'YScale','log')
    if r==3, xlabel('$\log_{10}(\varepsilon\, [\mathrm{W}/\mathrm{kg}])$',...
                    'Interpreter','latex'); end
    ylabel('PDF')
    yticks = -6:2:0;
    ytick_vals = 10.^yticks;
    set(gca, 'YTick', ytick_vals, ...
             'YTickLabel', arrayfun(@(x) sprintf('$10^{%d}$', x), yticks, 'UniformOutput', false), ...
             'TickLabelInterpreter','latex');
    set(gca,'Fontsize',fontsize,'fontname','times')
    box off;
end

% ---------- 3) row-level titles (annotation boxes) -----------------------
yPos = [0.93 0.63 0.34];        % vertical placements for three rows
for k = 1:3
   annotation('textbox',[0.13 yPos(k) 0.74 0.04],...
      'String',rowTitles{k},...
      'EdgeColor','none','HorizontalAlignment','center',...
      'FontSize',fontsize+2,'FontWeight','bold','FontName','times');
end

print(gcf,'-djpeg','-r200','Figure1.jpg')




% -------------------------------------------------------------------------
%  MLE local fitter + CDF helper (kept at end of file for clarity)
% -------------------------------------------------------------------------


