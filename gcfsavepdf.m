%% gcfsavepdf
% Save a figure to a PDF of appropriate size
%
% Tim Cronin // Massachusetts Institute of Technology
%
%%% Syntax
%   gcfsavepdf('/path/to/myfigure.pdf')
%
%%% Description
% Save a figure to a PDF. Before saving, this function queries the figure
% for its dimensions and set the PDF page to the appropriate size.
%
%%% Inputs
% *filename - PDF filename:*
% The filename (and, optionally, the path) where the PDF is saved.
%
function [ ] = gcfsavepdf( filename )

    set(gcf,'Units','Inches');
    figpos = get(gcf, 'Position');
    paperWidth = figpos(3);
    paperHeight = figpos(4);
    set(gcf, 'paperunits', 'Inches');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);

    print('-dpdf',filename);

end

