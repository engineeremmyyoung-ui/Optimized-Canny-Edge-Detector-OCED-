% THIS CODE SERVES DATA INTERACTIVE VISUALISATION TOOL FOR OCED PERFORMANCE
% EVALUATION METRICS ACROSS THE FOUR MEDICAL IMAGE MODALITIES

function OCED_App

% ---- Define metrics ----
metricNames = {'Accuracy','F1_score','PSNR','SSIM','Kappa','Precision','Recall','Entropy','Exec_time','MBE','MSE'};
modalities = {'X-ray','CT','MRI','Ultrasound'};

% Sample data
metricsData.Accuracy   = [0.9937, 0.9908, 0.9951, 0.9926];
metricsData.F1_score   = [0.9648, 0.9441, 0.9706, 0.9691];
metricsData.SSIM       = [0.9531, 0.9365, 0.9689, 0.9523];
metricsData.PSNR       = [22.981, 20.848, 24.404, 21.845];
metricsData.Kappa      = [0.9924, 0.9892, 0.9943, 0.9905];
metricsData.Precision  = [0.9679, 0.9690, 0.9769, 0.9702];
metricsData.Recall     = [0.9638, 0.9228, 0.9657, 0.9689];
metricsData.Entropy    = [0.4557, 0.4040, 0.4247, 0.5405];
metricsData.Exec_time  = [43.6, 34.9, 35.6, 41.8];
metricsData.MBE        = [0.0004, -0.0040, -0.0008, 0.0001];
metricsData.MSE        = [0.0063, 0.0092, 0.0049, 0.0074];

% ---- Create figure ----
f = figure('Name','OCED Metrics Viewer','Units','normalized','Position',[0.2 0.2 0.6 0.6]);

% ---- Metric Radiobuttons ----
bg = uibuttongroup('Parent',f,'Units','normalized','Position',[0.02 0.15 0.25 0.8]);
rbX = gobjects(length(metricNames),1);
for i = 1:length(metricNames)
    rbX(i) = uicontrol(bg,'Style','radiobutton','String',metricNames{i},...
        'Units','normalized','Position',[0.02, 1-(i*0.08), 0.95, 0.08],'HandleVisibility','off');
end

% Axes
ax = axes('Parent',f,'Units','normalized','Position',[0.32 0.15 0.65 0.75]);

% Attach radio button callback
bg.SelectionChangedFcn = @(src,evt) updatePlot(evt.NewValue.String, metricsData, ax, modalities, f);

% ---- Modality checkboxes ----
modalityCheckboxes = gobjects(length(modalities),1);
for i = 1:length(modalities)
    modalityCheckboxes(i) = uicontrol('Style','checkbox','String',modalities{i},...
        'Value',1,'Units','normalized',...
        'Position',[0.03+0.12*(i-1) 0.05 0.12 0.05],...
        'Callback',@(s,~) updatePlot(bg.SelectedObject.String, metricsData, ax, modalities, f));
end

% ---- Initial plot ----
updatePlot(metricNames{1}, metricsData, ax, modalities, f);

end

%% ---- Update Plot Function ----
function updatePlot(metricName, metricsData, ax, modalities, figHandle)
    % Collect which modalities are checked in the current figure only
    ch = findobj(figHandle,'Type','uicontrol','Style','checkbox');

    selected = [];
    for k = 1:length(ch)
        if ch(k).Value == 1
            idx = find(strcmp(modalities, ch(k).String));
            selected = [selected, idx];  % concatenate safely
        end
    end

    if isempty(selected)
        selected = 1:length(modalities); % default all
    end

    % Get metric values
    dataAll = metricsData.(metricName);
    y = dataAll(selected);

    % Clear axes and plot
    cla(ax)
    bar(ax, y, 'FaceColor',[0.2 0.6 0.8],'EdgeColor','k')
    set(ax,'XTick',1:length(y),'XTickLabel',modalities(selected),'FontSize',12)

    ylabel(ax,metricName,'FontSize',14)
    title(ax,['OCED ', metricName, ' per Modality'],'FontSize',16)
    grid(ax,'on')

    % Add numeric labels
    yl = ylim(ax);
    for ii = 1:length(y)
        text(ax,ii,y(ii)+0.02*(yl(2)-yl(1)),num2str(y(ii),'%.4g'),'HorizontalAlignment','center');
    end
end
