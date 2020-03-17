function [pks,loc] = f_peaks_PT(X,A,V,cps,win,plt,H,W)
% - find peaks or valleys in a profile transect
% 
% USAGE: [pks,loc] = f_peaks_PT(X,A,V,cps,win,plt,H,W)
%
% X   = structure of 'otolith profile' surface transect data created by f_cps2ppm_PT
% A   = cell array specifying the analyte to process;
%       ratios are used when 2 are provided e.g., A = {'Sr88' 'Ba137'};
% V   = find valleys instead of peaks                              (default = 0)
% cps = use raw cps data instead of ppm                            (default = 0)
% win = window size to filter/smooth time-series data             (default = 11)
% plt = create a plot                                              (default = 1)
% H   = height multplier                                           (default = 1)
% W   = width multplier                                            (default = 1)
% 
% pks = peaks or valleys
% loc = index to peaks or valleys

% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), V   = 0;  end % default don't find valleys
if (nargin < 4), cps = 0;  end % default don't use cps
if (nargin < 5), win = 11; end % default 11 point smoothing window
if (nargin < 6), plt = 1;  end % default create a plot
if (nargin < 7), H   = 1;  end % default height multplier is 1
if (nargin < 8), W   = 1;  end % default width multiplier is 1


% Check input:
if (numel(A)>2)
   error('A must specify only 1 or 2 analytes!')
end
if (numel(A)==2)
   ratio = 1;
else
   ratio = 0;
end
% -------------------------------------

% Extract fields from input:
iso = X.iso(:);
txt = X.txt(:);
pos = X.pos;
if (cps==1)
   Y = X.cps; % use raw cps data
else
   Y = X.ppm;
end

% Add combined element + isotope labels:
txt_iso = cellstr(strcat(txt,regexprep( cellstr(num2str(iso)),'\s', '') ))';

% Get index of analyte to process:
[TF,LOC] = ismember(A,txt_iso);
idx = find(TF==0);
if ~isempty(idx)
   fprintf('Target analytes were not found in A\n');
   fprintf('  %s \n\n',A{idx})
   error('All analytes in A must be found in X!');
end

% Retain only target analytes:
txt_iso = txt_iso(LOC);
Y       = Y(:,LOC);

% Optionally filter/smooth times-series data:
if (win>0)
   Y = f_filterSinclair(Y,1,win);
end

% Calculate ratios:
if (ratio == 1)
   Y    =  Y(:,1) ./ Y(:,2);
   yTxt = ['Ratio of ' txt_iso{1} '/' txt_iso{2}];
else
   yTxt = txt_iso;
end

% Set peak finding parameters:
H_var = 0.25*(max(Y) - min(Y)); % 25% of overall amplitude
W_var = ceil((win-1)/2);        % based on size of smoothing window

% Apply Height/Width multipliers:
H_var = H_var * H;
W_var = W_var * W;

% Set peak finding parameters:
% H_var = 0.03*(max(Y) - min(Y));
% W_var = 0.03*numel(Y);

if (V==1) % Find Valleys:
   [pks,loc] = findpeaks(max(Y) - Y,'MinPeakHeight',H_var,'MinPeakWidth',W_var);
   
   if (plt>0)
      figure;
      plot(pos,Y,pos(loc),max(Y)-pks,'or')
      title('Find Valleys')
      xlabel('Position (um)')
      ylabel(yTxt)
      % axis tight
   end
else % Find Peaks:
   [pks,loc] = findpeaks(Y,'MinPeakHeight',H_var,'MinPeakWidth',W_var);
   
   if (plt>0)
      figure;
      plot(pos,Y,pos(loc),pks,'or')
      title('Find Peaks')
      xlabel('Position (um)')
      ylabel(yTxt)
      % axis tight
   end
end
