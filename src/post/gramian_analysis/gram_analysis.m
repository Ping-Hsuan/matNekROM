% This script performs various plots related to Gramian matrix eigenvalues.
% Ensure that 'plotting_config.m' is in the same directory or provide the correct path.

close all; clear all;

% Load plotting configurations from 'plotting_config.m'
plotting_config; global myfig

%get_eig;

% Load eigenvalues from the 'eig.txt' file
eig = dlmread("eig.txt");

PlotGramEig(eig);           % Plot the Gramian matrix eigenvalue
PlotGramEigNormalized(eig); % Plot the normalized eigenvalue
PlotAccGramEig(eig, 0.9);   % Plot the accumulated eigenvalue
PlotEneCriteria(eig);       % Plot the energy critieria


function PlotEneCriteria(eig)
    % Plot the energy criteria: 1 - accumultaed eigenvalue
    % Inputs: 
    % Gramian matrix eigenvalue 
    global myfig

    f = figure(5);
    f = setFigureProperties(f, myfig); ax6 = axes(f);
    %
    loglog(ax6,1-(cumsum(eig)/sum(eig)),'-',myfig.lw,1); hold on
    yline(ax6,1e-1,'-'); hold on
    yline(ax6,1e-2,'-'); hold on
    yline(ax6,1e-3,'-'); hold on
    yline(ax6,1e-4,'-'); hold on
    %
    setAxisProperties(ax6, "$1-\sum^i_{j=1}\lambda_j/\sum^K_{j=1}\lambda_j$", myfig);
    print(gcf,"ene_criteria","-dpdf","-r300")
end

function PlotAccGramEig(eig, target_value)
    % Plot the accumulated eigenvalue of the gramian matrix
    % Inputs: 
    % Gramian matrix eigenvalue & target percentage
    global myfig

    f = figure(6);
    f = setFigureProperties(f, myfig); ax6 = axes(f);
    %
    acc_eig = cumsum(eig)/sum(eig);
    semilogx(ax6,acc_eig,myfig.lw,1.5); hold on;
    %
    xl = yline(ax6,target_value,'-',target_value*100+"% energy",'HandleVisibility','off'); hold on
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'right';
    xl.FontSize = 8;
    xl.LineWidth = 1.5;
    %
    % Find the index of the closest value
    [~, index] = min(abs(acc_eig-target_value));
    % Get the corresponding value
    closest_value = acc_eig(index);
    % Display the results
    fprintf('Index: %d\nValue: %f\n', index, closest_value);
    %
    plot(ax6,index, closest_value, 'ro', myfig.ms, 6);
    %
    setAxisProperties(ax6, "$\sum^i_{j=1}\lambda_j/\sum^K_{j=1}\lambda_j$", myfig);
    print(gcf,"accumulated_eig","-dpdf")
end


function fig = setFigureProperties(fig, myfig)
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'Units', 'Inches', 'Position', [0, 0, myfig.fig_width, myfig.fig_height], 'PaperUnits', 'Inches', 'PaperSize', [myfig.fig_width, myfig.fig_height]);
end

function setAxisProperties(ax, ylabelStr, myfig)
    ax.FontSize = 8;
    ax.YLabel.Interpreter = myfig.ltx;
    ax.YLabel.String = ylabelStr;
    ax.YLabel.FontSize = 8;
    ax.XLabel.Interpreter = myfig.ltx;
    ax.XLabel.String = "$i$";
    ax.XLabel.FontSize = 8;
    grid on;
    formatfig(ax);
end

function PlotGramEig(eig)
    % Plot the eigenvalue of the gramian matrix
    % Inputs: 
    % Gramian matrix eigenvalue 
    global myfig

    f  = figure(1); f  = setFigureProperties(f, myfig);
    ax = axes(f);
    %
    loglog(ax,eig,myfig.lw,1.5); hold on;
    %
    setAxisProperties(ax, "$\lambda_i$", myfig);
    print(gcf,"eig","-dpdf","-r300")
end

function PlotGramEigNormalized(eig)
    % Plot the normalized eigenvalue of the gramian matrix
    % Inputs: 
    % Gramian matrix eigenvalue 
    global myfig

    f  = figure(2);
    f  = setFigureProperties(f, myfig);
    ax = axes(f);
    %
    loglog(ax,eig/sum(eig),myfig.lw,1.5); hold on;
    %
    setAxisProperties(ax, "$\lambda_i/\sum^K_j \lambda_j$", myfig);
    %ylim([1e-12 1]);
    print(gcf,"norm_eig","-dpdf","-r300")
end