clear all; close all;

%% Load reduced mass matrix for computing L^2 norm
bu0 = dlmread("../../../ops/bu"); mb = sqrt(length(bu0)); bu0 = reshape(bu0,mb,mb);

%% Setting up plots and flags
plotting_config; setup_flag;
global iffom

%% Compute quantities for the ROM projection with N=mb
%% If one wants to compute at other N values, need to modify
%% snap -> snap(1:nb+1,:) and bu -> bu(1:nb+1,1:nb+1) 
snap = struct; 
snap.nb = mb-1;
snap.ucoef = dlmread("../../../ops/uk"); % Load projected coefficients of snapshots
ns = length(snap.ucoef)/mb; snap.ucoef = reshape(snap.ucoef,mb,ns);
snap.ua = mean(snap.ucoef,2); % Compute the mean of the projected coefficients
[snap.qoi, snap.mqoi, snap.stdqoi] = ConstructQoI(snap,bu0,ifene);

%% Reading FOM data if available
%% fomdir is defined in setup_flag
if iffom
    fom = struct();
    fom.data   = dlmread(fomdir);
    fom.t      = fom.data(:,2);
    fom.qoi    = fom.data(:,3);
    fom.mqoi   = mean(fom.qoi);
    fom.stdqoi = std(fom.qoi);
end

%% Plot quantities of the G-ROM
nb_list = [300];
grom = cell(1, length(nb_list));
if (ifgrom)
    for ii=1:size(nb_list,2)
        grom{ii} = SetupROMStructure();
        grom{ii}.nb = nb_list(ii);
    
        [grom{ii}.ucoef, grom{ii}.ua] = loadROMData("g-rom", grom{ii}.nb);
        [grom{ii}.qoi, grom{ii}.mqoi, grom{ii}.stdqoi] = ConstructQoI(grom{ii},bu0,ifene);
        grom{ii} = ComputeQoIError(grom{ii},snap,fom)
        DisplayResults(grom{ii},snap,fom,iffom);
        grom{ii} = GenerateTable(grom{ii},snap,fom,iffom);
    end
end

%% Plot quantities of the Regularized ROM
global ifreg
ifreg = 1;
reg_case = "EFR";
nb_list = [300];
regrom = cell(1, length(nb_list));
if (ifreg)
    [chi_list, radius_list] = create_hyperparam();
    for ii=1:size(nb_list,2)
        regrom{ii} = SetupROMStructure();
        regrom{ii}.nb = nb_list(ii);
        for order=4:4
            regrom{ii}.dforder=order;

            % Setup output directory
            outputdir = reg_case+"_N"+regrom{ii}.nb+"_"+regrom{ii}.dforder;
            mkdir(outputdir);
            for nn=1:size(chi_list,1)
                regrom{ii}.relax  = chi_list(nn);
                regrom{ii}.dfradius = radius_list(nn);
                disp([regrom{ii}.dforder regrom{ii}.dfradius regrom{ii}.relax])

                % Setup directory name
                casedir= sprintf('%s_%d_%d_%.8f_%.8f',reg_case,regrom{ii}.nb,regrom{ii}.dforder,regrom{ii}.dfradius,regrom{ii}.relax)
    
                [regrom{ii}.ucoef, regrom{ii}.ua] = loadROMData(casedir, regrom{ii}.nb);
                [regrom{ii}.qoi, regrom{ii}.mqoi, regrom{ii}.stdqoi] = ConstructQoI(regrom{ii},bu0,ifene);
                regrom{ii} = ComputeQoIError(regrom{ii},snap,fom)
                DisplayResults(regrom{ii},snap,fom,iffom);
                regrom{ii} = GenerateTable(regrom{ii},snap,fom,iffom);

                figure(1)
                set(gcf, 'PaperUnits', 'inches');
                set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
                    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height]);
                t=linspace(T_rom/size(regrom{ii}.qoi,1),T_rom,size(regrom{ii}.qoi,1))+T_0;
                plot(t,regrom{ii}.qoi,'-',cr,cmap(ii+1,:),dispname,getLegendLabel(reg_case,regrom{ii}.nb),lw,1.2); hold on
                plot(t,grom{ii}.qoi,'-',cr,cmap(1,:),dispname,"G-ROM with $N="+grom{ii}.nb+"$",lw,1.2); hold on
                xl = xline(T_snap,':',{'Training','window'},'HandleVisibility','off');
                xl.LabelVerticalAlignment = 'top';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize = 8;
%               xl.LabelOrientation = 'horizontal';
                xl.LineWidth = 1.5;

                if iffom
                    if if3dlidh && ifpred
                        t_fom=linspace(2725.125,3725,8000);
                        plot(t_fom,fom(:,3),'k-',lw,1.2); hold on
                    else
                        plot(fom.t,fom.qoi,'k-',lw,1.2,dispname,"FOM"); hold on
                    end
                else
                    t=linspace(T_snap/size(ene_snap,1),T_snap,size(ene_snap,1))+T_0;
                    plot(t,ene_snap,'k-',dispname,"Projection, $N="+nb+"$"); hold on
                end
                ax=gca; ax.FontSize=8; %xlim([0, T_rom])
                
                if ifene
                    if ~isempty(ylims)
                        ylim(ylims); 
                    end
                    if ~isempty(yticks_)
                        yticks(yticks_);
                    end
                    %ylim([1.1 1.15])
                    %ylim([3.1 3.3]);
                    %ylim([0.07 0.074]); yticks([0.07 0.071 0.072 0.073 0.074])
                    %ylim([0.015 0.0165]); %yticks([0.07 0.071 0.072 0.073 0.074])
                    xlabel("$t$",intp,ltx,fs,8); ylabel("$E(t)$",intp,ltx,fs,8);
                else
                    %ylim([1e-3, 1e-2]);
                    if ~isempty(ylims)
                        ylim(ylims); 
                    end
                    if ~isempty(yticks_)
                        yticks(yticks_);
                    end
                    %ylim([8e-4, 6e-3]); yticks([8e-4 2e-3 4e-3 6e-3])
                %   ylim([2e-4 2e-3]); yticks([2e-4 6e-4 1e-3 2e-3]) % 3dlid
                    %ylim([0.05 0.2]) % 3dminimal
                    xlabel("$t$",intp,ltx,fs,8); ylabel("$\mathrm{E}_{\mathrm{fluc}}(t)$",intp,ltx,fs,8);
                end
                %xlim([3000 4000])
                %if iffom
                %    LH(1) = plot(nan, nan, 'k-',lw,1.2);
                %    L{1} = "FOM";
                %else
                %    LH(1) = plot(nan, nan, 'k-',lw,1.2);
                %    L{1} = "Projection, $N=300$";
                %end
                %for ii=1:size(nb_list,2)
                %    nb=nb_list(ii);
                %%   LH(ii+1) = plot(nan, nan, '-.',cr,cmap(ii,:));
                %    LH(ii+1) = plot(nan, nan, '--',cr,cmap(ii,:),lw,1.2);
                %    L{ii+1} = "$N=$"+nb;
                %end
                
                figure(1)
                leg = legend({}, fs,8,intp,ltx,'location','best','NumColumns',2);
                leg.ItemTokenSize = [12,18];
                formatfig(ax); 
                
                casedir = strrep(casedir,".","p");
                if ifene
                    print(gcf,outputdir+"/"+"ene"+"_"+casedir,"-dpdf","-r300");
                else
                    print(gcf,outputdir+"/"+"intke"+"_"+casedir,"-dpdf","-r300");
                end
                close(1)
            end
            if (ifcreate_table)
                if ifene
                    if iffom
                        summary_reg.Properties.VariableNames = {'nb','order','radius','relax','fom.mqoi','fom.stdqoi','mene_snap','stdene_snap',...
                        'mene_rom','stdene_rom','mene_err','stdene_err','mene_err_f','stdene_err_f'};
                    else
                        summary_reg.Properties.VariableNames = {'nb','order','radius','relax','mene_snap','stdene_snap',...
                        'mene_rom','stdene_rom','mene_err','stdene_err'};
                    end
                    writetable(summary_reg,outputdir+"/"+"ene_grom"+"_"+reg_case+"_N"+reg.nb{ii}+"_"+dfOrder);
                else
                    if iffom
                        summary_reg.Properties.VariableNames = {'nb','order','radius','relax','mtke_fom','stdtke_fom','mtke_snap','stdtke_snap',...
                        'mtke_rom','stdtke_rom','mtke_err','stdtke_err','mtke_err_f','stdtke_err_f'};
                    else
                        summary_reg.Properties.VariableNames = {'nb','order','radius','relax','mtke_snap','stdtke_snap',...
                        'mtke_rom','stdtke_rom','mtke_err','stdtke_err'};
                    end
                    writetable(summary_reg,outputdir+"/"+"intke_grom"+"_"+reg_case+"_N"+reg.nb{ii}+"_"+dfOrder);
                end
            end
        end
    end
end

%% Plot quantities of the Constrained ROM
ifcrom = 0;
nb_list = [100];
crom = cell(1, length(nb_list));
if (ifcrom)
    for ii=1:size(nb_list,2)
        crom{ii} = SetupROMStructure();
        crom{ii}.nb = nb_list(ii);
        outputdir = "crom_N"+crom{ii}.nb;
        mkdir(outputdir);
    
        % Process C-ROM data
        [crom{ii}.ucoef, crom{ii}.ua] = loadROMData("c-rom", crom{ii}.nb);
        [crom{ii}.qoi, crom{ii}.mqoi, crom{ii}.stdqoi] = ConstructQoI(crom{ii},bu0,ifene);
        crom{ii} = ComputeQoIError(crom{ii},snap,fom)
        DisplayResults(crom{ii},snap,fom,iffom);
        crom{ii} = GenerateTable(crom{ii},snap,fom,iffom);

        figure(1)
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
            'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height]);
        t=linspace(T_rom/size(crom{ii}.qoi,1),T_rom,size(crom{ii}.qoi,1))+T_0;
        plot(t,crom{ii}.qoi,'-',cr,cmap(ii+1,:),dispname,"C-ROM with $N="+crom{ii}.nb+"$",lw,1.2); hold on
        if (ifgrom)
            plot(t,grom{ii}.qoi,'-',cr,cmap(1,:),dispname,"G-ROM with $N="+grom{1}.nb+"$",lw,1.2); hold on
        end
        xl = xline(T_snap,':',{'Training','window'},'HandleVisibility','off');
        xl.LabelVerticalAlignment = 'top';
        xl.LabelHorizontalAlignment = 'left';
%       xl.LabelOrientation = 'horizontal';
        xl.LineWidth = 1.5;

        if iffom
            if if3dlidh && ifpred
                t_fom=linspace(2725.125,3725,8000);
                plot(t_fom,fom(:,3),'k-',lw,1.2); hold on
            else
                plot(fom.t,fom.qoi,'k-',lw,1.2,dispname,"FOM"); hold on
            end
        else
            t=linspace(T_snap/size(ene_snap,1),T_snap,size(ene_snap,1))+T_0;
            plot(t,ene_snap,'k-',dispname,"Projection, $N="+nb+"$"); hold on
        end
        ax=gca; ax.FontSize=8; %xlim([0, T_rom])
        if ifene
            if ~isempty(ylims)
                ylim(ylims); 
            end
            if ~isempty(yticks_)
                yticks(yticks_);
            end
            xlabel("$t$",intp,ltx,fs,8); ylabel("$E(t)$",intp,ltx,fs,8);
        else
            if ~isempty(ylims)
                ylim(ylims); 
            end
            if ~isempty(yticks_)
                yticks(yticks_);
            end
            xlabel("$t$",intp,ltx,fs,8); ylabel("$\mathrm{E}_{\mathrm{fluc}}(t)$",intp,ltx,fs,8);
        end
        figure(1)
        leg = legend({}, fs,8,intp,ltx,'location','best','NumColumns',2);
        leg.ItemTokenSize = [12,18];
        formatfig(ax); 
        
        if ifene
            print(gcf,outputdir+"/"+"ene","-dpdf","-r300");
        else
            print(gcf,outputdir+"/"+"intke","-dpdf","-r300");
        end
        close(1)
    end
end

function [qoi, mqoi, stdqoi] = ConstructQoI(rom,bu0,ifene)
    bu = bu0(1:rom.nb+1,1:rom.nb+1);

    if ifene
        [qoi] = recon_ene(rom.ucoef,bu);
    else
        rom
        [qoi] = recon_intke(rom.ucoef,rom.ua,bu);
    end
    mqoi   = mean(qoi);
    stdqoi = std (qoi);
end

function [chi_list, radius_list] = create_hyperparam()
    chi_list = []; radius_list = [];
    chi = 0.005;
    chi_tmp = [chi];
    chi_tmp = [0.05 0.1 0.5 1];
    chi_tmp = [1];
    chi_tmp = [0.005];
    ttt = [linspace(0.1,0.3,20)'];
%   ttt = [0.13157895];
%   ttt = [0.25789474];
%   ttt = [0.16315789];
%   ttt = [linspace(0.01,0.1,25)';linspace(0.1,0.3,20)'];
    ttt = unique(ttt,'rows');
    for k=1:size(chi_tmp,2)
        chi = chi_tmp(k)
        for nn=1:size(ttt,1)
            chi_list= [chi_list;chi];
            radius_list = [radius_list;ttt(nn)];
        end
    end
end

function [ene] = recon_ene(ucoef,bu)
    % Recon energy from ucoef and bu
    ns = size(ucoef,2);
    prod_matrix = 0.5 * (ucoef.* (bu * ucoef));
    ene = sum(prod_matrix, 1)';
end

function label = getLegendLabel(reg_case, nb)
    if strcmp(reg_case, 'leray')
        label = "L-ROM with $N=" + nb + "$";
    else
        label = reg_case + "-ROM with $N=" + nb + "$";
    end
end

function [intke] = recon_intke(ucoef, ucoef_mean, bu)
    % Assuming ucoef is an nb x ns matrix
    % and bu is an nb x nb matrix

    % Calculate the term involving ucoef
    term_ucoef = sum(ucoef .* (bu * ucoef), 1)';

    % Calculate the term involving ucoef_mean
    term_ucoef_mean = 2 * sum(ucoef .* (bu * ucoef_mean), 1)';

    % Calculate the term involving ucoef_mean squared
    term_ucoef_mean_squared = sum(ucoef_mean .* (bu * ucoef_mean), 1)';
                                                                                                                                                                                                                                              
    % Combine the terms
    intke = term_ucoef - term_ucoef_mean + term_ucoef_mean_squared;
end

function [err] = computeRelativeError(ref, data)
    err = abs(ref - data) / ref;
end

function DisplayResults(rom,snap,fom,iffom)
    global ifreg

    if ifreg
        disp([rom.nb rom.dforder rom.dfradius rom.relax snap.mqoi snap.stdqoi rom.mqoi rom.stdqoi rom.merr rom.stderr]);
    else
        disp([rom.nb snap.mqoi snap.stdqoi rom.mqoi rom.stdqoi rom.merr rom.stderr]);
    end
    if iffom
        disp([rom.nb fom.mqoi fom.stdqoi rom.mqoi rom.stdqoi rom.merrF rom.stderrF]);
    end
end

function [rom] = SetupROMStructure();
    rom = struct();
    rom.ucoef = {};
    rom.ua = {};
    rom.nb = {};
    rom.qoi = {};
    rom.mqoi = {};
    rom.stdqoi = {};
    rom.merr = {};
    rom.stderr = {};
    rom.merrF = {};
    rom.stderrF = {};
    rom.summary = table;

    global ifreg
    if ifreg
        rom.dforder  = {};
        rom.dfradius = {};
        rom.relax    = {};
    end
end

function [ucoef, ua] = loadROMData(model, nb)
    global ifreg

    if ifreg
        dirname = "../" + model;
    else
        dirname = "../" + model + "_" + nb;
    end
    ucoef = dlmread(dirname + "/ucoef");
    ua = dlmread(dirname + "/ua");
    ndata = length(ucoef) / (nb + 1);
    ucoef = reshape(ucoef, ndata, nb + 1)';
end

function rom = GenerateTable(rom,snap,fom,iffom)
    global ifreg
    if ifreg
        results_cell = {rom.nb rom.dforder rom.dfradius rom.relax snap.mqoi snap.stdqoi rom.mqoi rom.stdqoi rom.merr rom.stderr};
    else
        results_cell = {rom.nb snap.mqoi snap.stdqoi rom.mqoi rom.stdqoi rom.merr rom.stderr};
    end

    if iffom
        results_cell = [results_cell, {fom.mqoi fom.stdqoi rom.merrF rom.stderrF}];
    end
    rom.summary = [rom.summary; results_cell];
end

function rom = ComputeQoIError(rom,snap,fom)
    global iffom
    [rom.merr]   = computeRelativeError(snap.mqoi, rom.mqoi);
    [rom.stderr] = computeRelativeError(snap.stdqoi, rom.stdqoi);

    if iffom
        [rom.merrF]   = computeRelativeError(fom.mqoi, rom.mqoi);
        [rom.stderrF] = computeRelativeError(fom.stdqoi, rom.stdqoi);
    end
end