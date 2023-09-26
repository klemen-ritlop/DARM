function plotDisplacements(ax, y, x, dy, dx, pointNames, selection)
% PLOTDISPLACEMENTS izriše geodetsko mrežo in premike točk na podane osi
%
% Inputs:
%    ax         : axes handle za koordinatne osi, na katerega želimo narediti izris
%    y          : vektor z y-koordinatami geodetskih točk
%    x          : vektor z x-koordinatami geodetskih točk
%    dx         : vektor komponent premikov geodetkih točk po x-osi
%    dy         : vektor komponent premikov geodetkih točk po y-osi
%    pointNames : vektor imen točk
%     selection : ali naj se posamezen premik izriše ali ne


    arguments
        ax {mustBeA(ax,'matlab.graphics.axis.Axes')}
        y {mustBeVector, mustBeNumeric, mustBeReal}
        x {mustBeVector, mustBeNumeric, mustBeReal}
        dy {mustBeVector, mustBeNumeric, mustBeReal}
        dx {mustBeVector, mustBeNumeric, mustBeReal}
        pointNames {mustBeVector, mustBeText}
        selection {mustBeVector, mustBeNumericOrLogical} = 1
    end

    cla(ax)

    if size(y, 2) ~= 1; y = y'; end
    if size(x, 2) ~= 1; x = x'; end
    if size(dy, 2) ~= 1; dy = dy'; end
    if size(dx, 2) ~= 1; dx = dx'; end
    if size(selection, 2) ~= 1; selection = selection'; end
    
    dy = dy .* selection;
    dx = dx .* selection;

    % meje mreže
    yMin = min(y);
    yMax = max(y);
    xMin = min(x);
    xMax = max(x);

    dyAbsMax = max(abs(dy));
    dxAbsMax = max(abs(dx));
    
    maxSize = max([yMax - yMin, xMax - xMin]);        % največja dimenzija mreže 
    m = (maxSize / 2) / max([dyAbsMax, dxAbsMax]);    % merilo vektorjev premikov
    
    % izračunaj koordinate koncev vektorjev deformacije (za izris v merilu m)
    y2 = y + dy * m;
    x2 = x + dx * m;                                            
    
    % meje mreža + vektorji 
    yMin = min([yMin, min(y2)]);
    yMax = max([yMax, max(y2)]);
    xMin = min([xMin, min(x2)]);
    xMax = max([xMax, max(x2)]);
    
    buffer = maxSize / 5;
    
    nPoints = length(y);
    
    % --- IZRIS ---
    hold(ax, 'on')

    % izriši mrežo
    for i = 1 : nPoints
        for j = i : nPoints
            plot(ax, [y(i), y(j)], [x(i), x(j)], '-', 'Color', '#e5e5e5');
        end
    end

    % izriši vektorje premikov
    for i = 1 : nPoints
        plot(ax, [y(i), y2(i)], [x(i), x2(i)], 'r-');
    end

    %izriši točke   
    plot(ax, y, x, 'ko');
    text(ax, y, x - buffer/7, pointNames, 'FontSize', 16)

    % nastavi osi
    ax.DataAspectRatio = [1, 1, 1];
    ax.XLim = [yMin - buffer, yMax + buffer];
    ax.YLim = [xMin - buffer, xMax + buffer];
    
    % izriši grafično merilo mreže
    e = floor(log10(maxSize/2));
    scalebarLength = 10^e;
    plot(ax, [yMin - buffer/2, yMin - buffer/2 + scalebarLength], [xMin - 3*buffer/4, xMin - 3*buffer/4] , '|-', 'Color', '#a0a0a0')
    text(ax, yMin - buffer/2 + scalebarLength + maxSize/50, xMin - 3*buffer/4 , [num2str(scalebarLength) ' m']);

    % izriši grafično merilo premikov
    e = floor(log10(max([dyAbsMax, dxAbsMax])));
    scalebarLength = 10^e;
    plot(ax, [yMin - buffer/2, yMin - buffer/2 + scalebarLength * m], [xMin - buffer/2, xMin - buffer/2] , 'r|-')
    text(ax, yMin - buffer/2 + scalebarLength * m + maxSize/50, xMin - 2*buffer/4 , [num2str(scalebarLength*1000) ' mm']);

    %l = line([text_end_h, text_end_h + m], [text_mid_v, text_mid_v]);
    %ax.Toolbar.Visible = 'off';
    disableDefaultInteractivity(ax)

    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    ax.XTickLabelRotation = 0;

end
