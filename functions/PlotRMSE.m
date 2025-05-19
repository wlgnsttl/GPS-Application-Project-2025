function PlotRMSE(TTs, NEV, vis_sats, used_sats)
    %   in  :   TTs:    gs time tag
    %           NEV:    dNEV 
    %           visible sat by epoch
    %           used sat by epoch
    %   out :   fig

    %% ── 1. 레이아웃 파라미터 블록 (여기만 수정) ────────────────────────
    % [가로 배치]
    leftX    = 0.06;    % 왼쪽 열 시작
    leftW    = 0.28;    % 왼쪽 열 폭
    hGap     = 0.05;    % 열 간 가로 간격
    midW     = 0.20;    % 가운데 열 폭
    midX     = leftX + leftW + hGap;          % 가운데 열 시작
    rightW   = 0.31;    % 오른쪽 열 폭 (이 값만 바꾸세요)
    rightX   = midX + midW + hGap;            % 오른쪽 열 시작

    % [세로 배치]
    botY      = 0.08;   % 바닥 여유
    topMargin = 0.08;   % 상단 여유
    contentH  = 1 - botY - topMargin;  % 콘텐츠 전체 높이
    vGap      = 0.10;   % 중간 그래프 세로 간격
    midH      = (contentH - vGap)/2;   % 가운데 두 그래프 각각 높이
    tableH    = 0.12;   % 오른쪽 하단 테이블 높이
    graphH    = contentH - tableH - vGap;  % 오른쪽 상단 그래프 높이

    % [테이블 내부]
    rowGap    = 0.25;   % 테이블 헤더 높이 비율
    fontSize  = 11;     % 테이블 글씨 크기
    % ────────────────────────────────────────────────────────────────

    %% 2. RMSE 계산

    [rmse, horizontal, vertical, rmse3d] = nev2rmse(NEV);
    [~, t] = gs2gdhms(TTs);

    rmse_v  = rmse(2);
    rmse_h  = rmse(1);
    rmse_3d = rmse(3);

    north = NEV(:,1);
    east = NEV(:,2);

    %% 3. Figure 생성
    figure;

    %% 4. 왼쪽 산점도
    axes('Position',[leftX, botY, leftW, contentH]);
    hold on; grid on; axis equal;
    scatter(east, north, 5, 'r','filled');
    theta = linspace(0,2*pi,200);
    round = max(abs(horizontal));
    plot((round/2)*cos(theta),(round/2)*sin(theta),'b','LineWidth',1.2);
    plot(round*cos(theta),round*sin(theta),'b','LineWidth',1.2);
    xlabel('East(m)'); ylabel('North(m)');
    title(sprintf('Circle Rad : %.2fm, %.2fm',round/2,round));
    xline(0,'LineWidth',1); yline(0,'LineWidth',1);
    xlim([-(round*1.2) (round*1.2)]); ylim([-(round*1.2) (round*1.2)]);

    %% 5. 가운데 상단: Vertical 오차
    axes('Position',[midX, botY+midH+vGap, midW, midH]);
    plot(t, vertical, 'k'); grid on;
    xlabel('Hours'); ylabel('Vertical (m)');
    title(sprintf('RMSE : %.4fm', rmse_v));
    % ylim([-10 10]);

    %% 6. 가운데 하단: Horizontal 오차
    axes('Position',[midX, botY, midW, midH]);
    plot(t, horizontal, 'k'); grid on;
    xlabel('Hours'); ylabel('Horizontal (m)');
    title(sprintf('RMSE : %.4fm', rmse_h));
    % ylim([0 7]);

    %% 7. 오른쪽 상단: 3D RMSE + 위성 수
    axes('Position',[rightX, botY+tableH+vGap, rightW, graphH]);
    yyaxis left;
        plot(t, rmse3d, 'k');
        ylabel('3D(m)'); 
        % ylim([0 14]);
    yyaxis right;
        stairs(t, vis_sats, 'm','LineWidth',1); hold on;
        stairs(t, used_sats,'b--','LineWidth',1);
        ylabel('Satellites'); ylim([0 15]);
    grid on;
    xlabel('Hours');
    legend({'rmse','Visible Sats','Used Sats'}, 'Location','northwest');
    title(sprintf('RMSE : %.4fm', rmse_3d));

    %% 8. 오른쪽 하단: 3열 테이블 (헤더 하단만 두껍게)
    axT = axes('Position',[rightX, botY, rightW, tableH], 'Visible','off');
    hold(axT,'on'); xlim(axT,[0 1]); ylim(axT,[0 1]);

    % 테두리
    plot(axT,[0 1],[1 1],'k','LineWidth',1);
    plot(axT,[0 1],[0 0],'k','LineWidth',1);
    plot(axT,[0 1],[1-rowGap 1-rowGap],'k','LineWidth',2);

    % 헤더 텍스트
    headers = {'V RMSE','2D(H) RMSE','3D RMSE'};
    for i = 1:3
        x = (i-0.5)/3;
        text(axT, x, 1-rowGap/2, headers{i}, ...
             'FontWeight','bold','FontSize',fontSize, ...
             'HorizontalAlignment','center');
    end

    % 데이터 텍스트
    yData = (1-rowGap)/2;
    values = {sprintf('%.2fm',rmse_v), sprintf('%.2fm',rmse_h), sprintf('%.2fm',rmse_3d)};
    for i = 1:3
        x = (i-0.5)/3;
        text(axT, x, yData, values{i}, ...
             'FontSize',fontSize, 'HorizontalAlignment','center');
    end

    axis(axT,'off');
end
