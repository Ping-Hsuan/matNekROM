close all; clear all;

hdr;

%get_eig;

eig = dlmread("eig.txt");

% Plot the eigenvalue of the gramian matrix
f=figure(1);
f = setFigureProperties(f, myfig); ax = axes(f);
%
loglog(ax,eig,lw,1.5); hold on;
%
setAxisProperties(ax, "$\lambda_i$", myfig);
print(gcf,"eig","-dpdf","-r300")

% Plot the normalized eigenvalue of the gramian matrix
f=figure(2);
f = setFigureProperties(f, myfig); ax = axes(f);
%
loglog(ax,eig/sum(eig),lw,1.5); hold on;
%
setAxisProperties(ax, "$\lambda_i/\sum^K_j \lambda_j$", myfig);
%ylim([1e-12 1]);
print(gcf,"norm_eig","-dpdf","-r300")

f=figure(5);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
loglog(1-(cumsum(eig)/sum(eig)),'-',lw,1); hold on
loglog([1 500], [1e-1 1e-1]); hold on
loglog([1 500], [1e-2 1e-2]); hold on
loglog([1 500], [1e-3 1e-3]); hold on
loglog([1 500], [1e-4 1e-4]); hold on
ax=gca;
ax.FontSize=6;
leg.ItemTokenSize = [10,18]
ylabel("$1-\sum^i_{j=1}\lambda_j/\sum^K_{j=1}\lambda_j$",intp,ltx,fs,8);
xlabel("$i$",intp,ltx,fs,8);
grid on
print(gcf,"ene_criteria","-dpdf","-r300")

% Plot the accumulated eigenvalue of the gramian matrix
f = figure(6);
f = setFigureProperties(f, myfig); ax6 = axes(f);
%
acc_eig = cumsum(eig)/sum(eig);
semilogx(ax6,acc_eig,lw,1.5); hold on;
%
target_value = 0.9;
xl = yline(ax6,target_value,'-','90% energy','HandleVisibility','off'); hold on
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize = 8;
xl.LineWidth = 1.5;
%
% Find the index of the closest value
[~, index] = min(abs(acc_eig- target_value));
% Get the corresponding value
closest_value = acc_eig(index);
% Display the results
fprintf('Index: %d\nValue: %f\n', index, closest_value);
%
plot(ax6,index, closest_value, 'ro', myfig.ms, 6);
%
setAxisProperties(ax6, "$\sum^i_{j=1}\lambda_j/\sum^K_{j=1}\lambda_j$", myfig);
print(gcf,"accumulated_eig","-dpdf")


function fig = setFigureProperties(fig, myfig)
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'Units', 'Inches', 'Position', [0, 0, myfig.fig_width, myfig.fig_height], 'PaperUnits', 'Inches', 'PaperSize', [myfig.fig_width, myfig.fig_height])
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