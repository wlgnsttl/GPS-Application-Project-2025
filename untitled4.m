fig = gcf;  % 또는 openfig('filename.fig','invisible');
axesList = findobj(fig, 'Type', 'axes');

% 결과 저장 구조
recovered = struct();

for i = 1:length(axesList)
    ax = axesList(i);
    lines = findobj(ax, '-property', 'YData');

    % 왼쪽 산점도 (axis equal + scatter + 제목에 'Circle')
    if strcmp(get(ax, 'DataAspectRatioMode'), 'manual') && contains(get(get(ax, 'Title'), 'String'), 'Circle')
        sc = findobj(ax, 'Type', 'scatter');
        recovered.NEV = [get(sc, 'YData')', get(sc, 'XData')'];  % North, East

    % 가운데 상단 (Vertical)
    elseif contains(get(get(ax, 'YLabel'), 'String'), 'Vertical')
        recovered.vertical = get(lines, 'YData')';
        recovered.t = get(lines, 'XData')';

    % 가운데 하단 (Horizontal)
    elseif contains(get(get(ax, 'YLabel'), 'String'), 'Horizontal')
        recovered.horizontal = get(lines, 'YData')';

    % 오른쪽 상단 (3D + Sats)
    elseif contains(get(get(ax, 'YLabel'), 'String'), '3D') || contains(get(get(ax, 'Title'), 'String'), '3D')
        allLines = findobj(ax, 'Type', 'line');
        for j = 1:length(allLines)
            y = get(allLines(j), 'YData')';
            x = get(allLines(j), 'XData')';
            style = get(allLines(j), 'LineStyle');
            color = get(allLines(j), 'Color');

            if strcmp(style, '-') && isequal(color, [0 0 0])  % 검정 실선 → rmse3d
                recovered.rmse3d = y;
            elseif strcmp(style, '--')  % 점선 → used_sats
                recovered.used_sats = y;
            elseif strcmp(style, '-') && isequal(color, [1 0 1])  % 자주색 실선 → vis_sats
                recovered.vis_sats = y;
            end
        end
    end
end
