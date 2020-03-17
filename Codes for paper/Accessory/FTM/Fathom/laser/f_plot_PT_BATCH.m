function f_plot_PT_BATCH(fname,A,type,sm,verb)
% - batch plot of profile transect data
%
% USAGE: f_plot_PT_RATIO(fname,A,type,sm,verb);
%
% fname = cell array of file(s) to process
% A   = cell array specifying 2 analytes to plot as a ratio
%       e.g., A = {'Sr88' 'Ba137'};
% type  = 1 = ppm, 2 = cps, 3 = ratio
% sm    = filter/smooth time-series data                           (default = 1)
% verb  = verbose output of progress                               (default = 1)

% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 4), sm    = 1; end % default filter/smooth times-series data
if (nargin < 5), verb  = 1; end % default verbose output of progress

% Check input:
if type==2
   if numel(A) ~=2
      error('For ratios A must specify 2 analytes!')
   end
end

% Specify Regular Expression matching pattern:
ptn  = '^tra_\d+$'; % trasect naming convention
% -------------------------------------

% Check that FNAME is a cell array:
if ~iscell(fname), error('FNAME must be a cell array'); end

% Check that files in FNAME exist and contain the proper variables:
nFiles = numel(fname); % get # files
for i = 1:nFiles
   % Check that file exists:
   if (exist(fname{i},'file')~=2)
      error(['Cannot find file ' fname{i} '!'])
   end
   
   % Check for presence of profile transects:
   if isempty(who('-file',fname{i},'-regexp',ptn));
      error(['File ' fname{i} ' contains no profile transects!']);
   end
end
% -------------------------------------

% Process each file separately:
for i = 1:nFiles
   
   % Show file being processed:
   if (verb>0)
      fprintf('\nProcessing file %s...\n',fname{i});
   end
   
   % Load variables:
   load(fname{i},'-regexp',ptn); % profile transects
   
   % Get list of transects:
   uT = who('-file',fname{i},'-regexp',ptn);
   nT = numel(uT);     % # transects in this file
   
   % Plot each transect separately:
   for j = 1:nT
      
      % Extract fields from input:
      X   = eval(uT{j});
      iso = X.iso(:);
      txt = X.txt(:);
      pos = X.pos;
      if type==2
         ppm = X.cps; % use raw cps data
      else
         ppm = X.ppm;
      end
      
      % Add combined element + isotope labels:
      txt_iso = cellstr(strcat(txt,regexprep( cellstr(num2str(iso)),'\s', '') ))';
      
      % Get index of analytes to plot:
      [TF,LOC] = ismember(A,txt_iso);
      idx = find(TF==0);
      if ~isempty(idx)
         fprintf('Target analytes were not found in %s:\n',fname{i});
         fprintf('  %s \n\n',A{idx})
         error('All analytes in A must be found in X!');
      end
      
      % Retain only target analytes:
      txt_iso = txt_iso(LOC);
      ppm     = ppm(:,LOC);
      
      % Optionally filter/smooth times-series data:
      if (sm>0)
         ppm = f_filterSinclair(ppm);
      end
      
      if (type<3) % PPM or CPS
         
         % Create plot:
         figure;
         set(gcf,'color','w'); % set background color to white
         hold on;
                  
         % Plot each analyte separately:
         nc   = size(ppm,2); % # analytes to plot
         hdl = NaN(1,nc);   % preallocate
         
         for k = 1:nc
            if nc>10
               hdl(k) = plot(pos,ppm(:,k),'Color',f_rgb(k),'LineStyle',f_style(l));
            else
               hdl(k) = plot(pos,ppm(:,k),'Color',f_rgb(k),'LineStyle','-');
            end
         end
         
         % Customize plot:
         title(uT{j}, 'Interpreter','none')
         xlabel('Position (um)')
         if type==1
            yTxt = 'Concentration (ppm)';
         else
            yTxt = 'raw CPS';
         end
         ylabel(yTxt);
         box on;
         
         % Create legend:
         legend(hdl,txt_iso);
         
      else % RATIOS
         % Plot analyte ratios:
         figure
         plot(pos,ppm(:,1)./ppm(:,2),'Color','k','LineStyle','-');
         
         % Customize plot:
         title(uT{j}, 'Interpreter','none')
         xlabel('Position (um)')
         yTxt = ['Ratio of ' txt_iso{1} '/' txt_iso{2}];
         ylabel(yTxt);
         box on;
      end
      
      % Save as PDF's:
      f_pdf(uT{j})
      
   end
   % Clear variables associated with this file:
   clear uT;
end
