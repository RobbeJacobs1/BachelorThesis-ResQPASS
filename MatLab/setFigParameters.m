function setFigParameters(filename, width, height, fig, padding)
% This script is used to set the size of all images for the paper.
    if nargin == 0
        fig = gcf;
        width = 15.5;
        height = 7.75;
        filename = [];
        padding = 10;
    elseif nargin == 1
        fig = gcf;
        width = 15.5;
        height = 7.75;
        padding = 10;
    elseif nargin == 3
        fig = gcf;
        padding = 10;
    elseif nargin == 4
        padding = 10;
    end
    if isempty(fig)
        fig = gcf;
    end
    if isempty(width)
        width = 15.5;
    end
    if isempty(height)
        height = 7.75;
    end

    ax = findall(fig,'type','axes');
    for i = 1:length(ax)
        set(ax(i), 'LooseInset', get(ax(i), 'TightInset'))
    end

    % Onscreen
    set(fig, "Units", "centimeters");
    set(fig, "Position", [0 0 width height]);
    
    % Export
    set(fig, "PaperUnits", "Centimeters");
    set(fig, "PaperSize", [width height]);
    set(fig, "PaperPosition", [0 0 width height])

    if ~isempty(filename)
        exportgraphics(fig, "Images/" + filename + ".pdf", ...
            "ContentType", "vector", ...
            "Padding", padding);
    end
end