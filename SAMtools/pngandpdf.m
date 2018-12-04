%% png2mp4
%
% Output PNG and PDF versions of a figure.
% 
% Tristan Abbott // MIT EAPS // 12/4/2018
function png2mp4(fname)

    print([fname '.png'], '-dpng', '-r200');
    gcfsavepdf([fname '.pdf']);

end
