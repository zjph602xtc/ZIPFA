function f_export_ppm(X,fname)
% - export LA-ICP-MS parsed spot PPM & LOD data
%
% USAGE: f_export_ppm(X,'fname')
%
% X     = structure of processed spot data by f_cps2ppm_SPOT
% fname = name of destination file (omit *.csv extension)
%
% SEE ALSO: f_export_PT, f_export_cps

% -----Author:-----
% by David L. Jones, Apr-2017
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
% Don't overwrite existing file:
if exist([fname '_ppm.csv'],'file')
   error('Destination file already exists!')
end
% -------------------------------------

% Export PPM data as CSV:
f_exportR([X.ppm],cellstr(num2str(X.oto)),...
   [X.txt_iso],[fname '_ppm.csv']);

% Export LOD data as CSV:
f_exportR([X.LOD],cellstr(num2str(X.oto)),...
   [X.txt_iso],[fname '_LOD.csv']);

