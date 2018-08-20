%% Generate figures with standard screen-independent sizes
function fig = standard_figure(ftype)

    fig = figure('Units', 'Inches', 'PaperUnits', 'Inches');
    switch ftype
        case 'medium'
            set(fig, 'Position', [1 1 8 6], 'PaperSize', [8 6]);
        case 'medium-wide'
            set(fig, 'Position', [1 1 12 6], 'PaperSize', [12 6]);
        case 'medium-narrow'
            set(fig, 'Position', [1 1 6 6], 'PaperSize', [6 6]);
        case 'wide'
            set(fig, 'Position', [1 1 12 6], 'PaperSize', [12 6]);
        case 'half-page'
            set(fig, 'Position', [1 1 3.75 2.82], 'PaperSize', [3.75 2.82]);
        case 'half-page-short'
            set(fig, 'Position', [1 1 3.75 1.5], 'PaperSize', [3.75 1.5]);
        otherwise
            error('Unrecognized standard figure: %s', ftype);
    end

end
