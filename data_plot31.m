%% Plot the reconfiguration ratio varying with cost coefficient
%%
line_width = [1 1 1 1 1];
color_set = [Color.Black; Color.MildGreen; Color.Purple; Color.MildBlue; Color.Red; ...
    Color.Black; Color.Gray];
legend_label = {'Baseline', 'RSCV', 'RSEV', 'HRSCV', 'HRSEV', ...
    'RSCV-FI', 'RSEV-FI', 'HRSCV-FI', 'HRSEV-FI'};
line_style = {':', '-.', 'none', '-', 'none', '--', 'none', ':', 'none'};
marker = {'none', 'none', 's', 'none', 'd', 'none', '+', 'none', 'x'};
marker_size = [6 6 6 6 6 6 7 6 7];

%% normalized
if exist('fig_stat_reconfig', 'var') && fig_stat_reconfig.isvalid
    figure(fig_stat_reconfig);
else
    fig_stat_reconfig = figure('Name', 'Satistics of Reconfiguration');
end
x = 1:length(etas);
ratio_reconfig = results.Reconfig{idx,'ReVariables'}./...
    results.Reconfig{idx,'NumVariables'};
avgrat_reconfig = mean(ratio_reconfig)*ones(length(x),1);
avgfair_reconfig = mean(results.Reconfig{idx,'FairIndex'})*ones(length(x),1);
avgrat_fastconfig = zeros(length(x),1);
avgrat_fastconfig2 = zeros(length(x),1);
avgfair_fastconfig = zeros(length(x),1);
avgfair_fastconfig2 = zeros(length(x),1);
if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
    avgrat_dimconfig = zeros(length(x),1);
    avgrat_dimconfig2 = zeros(length(x),1);
    avgfair_dimconfig = zeros(length(x),1);
    avgfair_dimconfig2 = zeros(length(x),1);
end
sel_index = 51:100;
for i = x
    ratio_fastconfig = results.Fastconfig{i}{idx,'ReVariables'}./...
        results.Fastconfig{i}{idx,'NumVariables'};
    ratio_fastconfig2 = results.Fastconfig2{i}{idx,'ReVariables'}./...
        results.Fastconfig2{i}{idx,'NumVariables'};
    avgrat_fastconfig(i) = mean(ratio_fastconfig);
    avgrat_fastconfig2(i) = mean(ratio_fastconfig2);
    avgfair_fastconfig(i) = mean(results.Fastconfig{i}{idx(sel_index),'FairIndex'});
    avgfair_fastconfig2(i) = mean(results.Fastconfig2{i}{idx(sel_index),'FairIndex'});
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        ratio_dimconfig = results.Dimconfig{i}{idx,'ReVariables'}./...
            results.Fastconfig{i}{idx,'NumVariables'};
        ratio_dimconfig2 = results.Dimconfig2{i}{idx,'ReVariables'}./...
            results.Dimconfig2{i}{idx,'NumVariables'};
        avgrat_dimconfig(i) = mean(ratio_dimconfig);
        avgrat_dimconfig2(i) = mean(ratio_dimconfig2);
        avgfair_dimconfig(i) = mean(results.Dimconfig{i}{idx(sel_index),'FairIndex'});
        avgfair_dimconfig2(i) = mean(results.Dimconfig2{i}{idx(sel_index),'FairIndex'});
    end
end
yyaxis left;
hl = plot(...
    etas, avgrat_fastconfig./avgrat_reconfig, '-.^',...
    etas, avgrat_fastconfig2./avgrat_reconfig, '-.s');
% hl = plot(etas, avgrat_reconfig, '--',...
%     etas, avgrat_fastconfig, '-.^',...
%     etas, avgrat_fastconfig2, '-.s');
if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
    hold on;
    hl = plot(etas, avgrat_dimconfig./avgrat_reconfig, '-x',...
        etas, avgrat_dimconfig2./avgrat_reconfig, '-+');
    hold off;
    hl = hl(1).Parent.Children;
end
nl = length(hl);
for k=1:nl   % REVERSE order
    hl(k).Color = color_set(nl-k+2).RGB;
    hl(k).LineWidth = line_width(nl-k+2);
    hl(k).Marker = marker{nl-k+2};
    hl(k).MarkerSize = marker_size(nl-k+2);
    hl(k).LineStyle = line_style{nl-k+2};
end
if length(etas)>=10
    marker_index = round(linspace(1,length(etas), 10));
    for k=1:nl
        hl(k).MarkerIndices = marker_index;
    end
end
xlim([etas(1), etas(end)]);
hl(1).Parent.YLim(1) = 0;
hl(1).Parent.YLim(2) = min(hl(1).Parent.YLim(2)+0.1, 1);
xlabel('Parameter \eta');
ylabel('Normalized reconfigration ratio');
yyaxis right
hr = plot(...
    etas, avgfair_fastconfig, '--^',...
    etas, avgfair_fastconfig2, '--s');
if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
    hold on;
    hr = plot(etas, avgfair_dimconfig, '--x',...
        etas, avgfair_dimconfig2, '--+');
    hold off;
    hr = hr(1).Parent.Children;
end
nr = length(hr);
for k=1:nl   % REVERSE order
    hr(k).Color = color_set(nl-k+2).RGB;
    hr(k).LineWidth = line_width(nl-k+2);
    hr(k).Marker = marker{5+k};
    hr(k).MarkerSize = marker_size(5+k);
    hr(k).LineStyle = line_style{5+k};
end
if length(etas)>=10
    marker_index = round(linspace(1,length(etas), 10));
    for k=1:nl
        hr(k).MarkerIndices = marker_index;
    end
end
legend(legend_label((1:(length(hl)+length(hr)))+1), 'Location', 'best', 'Orientation', 'vertical');
ylabel('Average fairness index');

%% original
if exist('fig_stat_num_reconfig', 'var') && fig_stat_num_reconfig.isvalid
    figure(fig_stat_num_reconfig);
else
    fig_stat_num_reconfig = figure('Name', 'Satistics of Reconfiguration');
end
x = 1:length(etas);
ratio_reconfig = results.Reconfig{idx,'ReVariables'}./...
    results.Reconfig{idx,'NumVariables'};
avgrat_reconfig = mean(ratio_reconfig)*ones(length(x),1);
avgrat_fastconfig = zeros(length(x),1);
avgrat_fastconfig2 = zeros(length(x),1);
if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
    avgrat_dimconfig = zeros(length(x),1);
    avgrat_dimconfig2 = zeros(length(x),1);
end
for i = x
    ratio_fastconfig = results.Fastconfig{i}{idx,'ReVariables'}./...
        results.Fastconfig{i}{idx,'NumVariables'};
    ratio_fastconfig2 = results.Fastconfig2{i}{idx,'ReVariables'}./...
        results.Fastconfig2{i}{idx,'NumVariables'};
    avgrat_fastconfig(i) = mean(ratio_fastconfig);
    avgrat_fastconfig2(i) = mean(ratio_fastconfig2);
    if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
        ratio_dimconfig = results.Dimconfig{i}{idx,'ReVariables'}./...
            results.Fastconfig{i}{idx,'NumVariables'};
        ratio_dimconfig2 = results.Dimconfig2{i}{idx,'ReVariables'}./...
            results.Dimconfig2{i}{idx,'NumVariables'};
        avgrat_dimconfig(i) = mean(ratio_dimconfig);
        avgrat_dimconfig2(i) = mean(ratio_dimconfig2);
    end
end
hl = plot(etas, avgrat_reconfig, '--',...
    etas, avgrat_fastconfig, '-.^',...
    etas, avgrat_fastconfig2, '-.s');
if ~exist('b_disable_dimconfig', 'var') || b_disable_dimconfig == false
    hold on;
    hl = plot(etas, avgrat_dimconfig, '-x',...
        etas, avgrat_dimconfig2, '-+');
    hold off;
end
hl = hl(1).Parent.Children;
nl = length(hl);
for k=1:nl   % REVERSE order
    hl(k).Color = color_set(nl-k+1).RGB;
    hl(k).LineWidth = line_width(nl-k+1);
end
if length(etas)>=10
    marker_index = round(linspace(1,length(etas), 10));
    for k=1:nl
        hl(k).MarkerIndices = marker_index;
    end
end
xlim([etas(1), etas(end)]);
hl(1).Parent.YLim(1) = 0;
hl(1).Parent.YLim(2) = min(hl(1).Parent.YLim(2)+0.1, 1);
legend(legend_label(1:length(hl)), 'Location', 'best', 'Orientation', 'vertical');
xlabel('Parameter \eta/T');
ylabel('Avrage reconfigration ratio');

%%
% If |MarkerIndices| is not availble for old version MATLAB, use the following code
% instead.
%
   MarkerIndices = [1, 5, 10];
   myplot = plot(x, y, 'b-.o');
   hold on;
   mymarkers = plot(x(MarkerIndices), y(MarkerIndices), 'ro');
   legend(myplot)