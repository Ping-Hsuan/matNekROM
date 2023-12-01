lw='linewidth';               %%% Plotting defs
fs='fontsize';                %%% Plotting defs
intp = 'interpreter';         %%% Plotting defs
ltx  = 'latex';               %%% Plotting defs
cr = 'Color';
ms = 'MarkerSize';
dispname = 'DisplayName';
colour_shape = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840]};

%%columnwidth = 469.75499; % in pts
columnwidth = 469.75499/2;
fig_width_pt = columnwidth;
inches_per_pt = 1.0/72.27;
fig_width = fig_width_pt*inches_per_pt;
heightratio = (sqrt(5)-1.0)/2.0;
%heightratio = 0.50;
heightratio = 0.40;
fig_height = fig_width*heightratio;
%fig_height_pt = (259.04362);
%fig_height = fig_height_pt*inches_per_pt;
fig_size = [fig_width, fig_height];

markers = ['o','x','square'];
linest = ["-","--","-."];

