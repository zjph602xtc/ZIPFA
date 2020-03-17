% Example of performing Variation Partitioning Analysis using 3 sets of
% explanatory matrices.
% 
% by David L. Jones, Apr-2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example follows the analysis presented in Legendre & Legendre (2012:
% p.658).
% 
% The file 'Doubs.xls' includes data collected by Verneaux (1973) at 30 sampling
% sites from the Doubs River near the border of Switzerland and France. The
% 'bio' tab contains the abundances of 27 species of fishes collected in the
% river; the 'env' tab contains 11 environmental variables related to the
% hydrology, geomorphology, and chemistry of the river; and the 'spatial' tab
% includes geographical coordinates (as Cartesian X & Y)  of the sampling
% sites. Abbreviations for the environmenal data are coded as follows;
% 
% das = distance from shore (km)
% alt = altitude (m)
% pen = slope (per mille)
% deb = mean minimum discharge (m^3 s^-1)
% ph  = pH of water
% dur = Calcium concentration (hardness) (mg L^-1)
% pho = Phosphate concentration (mg L^-1)
% nit = Nitrate concentration ((mg L^-1)
% amm = Ammonium concentration (mg L^-1)
% oxy = Dissolved oxygen (mg L^-1)
% dbo = Biological oxygen demand (mg L^-1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace:
clz;

% Import the data:
bio = f_importCSV('Doubs_bio.csv',1);     % fish abundances
env = f_importCSV('Doubs_env.csv',1);     % environmental data

% Check the biological data to see if any sites have no fish present:
find(sum(bio.dat,2)==0)
% 
% ans =
% 
%      8

% Remove site 8 from the data set:
bio.dat(8,:) = [];
bio.txt(8)   = [];
env.dat(8,:) = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partition the environmental data from the Doubs River into three sets of
% variables: Topography (alt, pen, and deb), Chemistry (ph, dur, pho, nit, amm,
% oxy, and dbo), and Geography (das).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create three subsets of the env data:
top.dat = env.dat(:,2:4);  % topography
top.txt = env.txt(2:4);

che.dat = env.dat(:,5:11); % chemistry
che.txt = env.txt(5:11);

geo.dat = env.dat(:,1);    % geography
geo.txt = env.txt(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Traditional RDA:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hellinger-transform the biotic data:
vpa = f_vpa(f_hellinger(bio.dat),top.dat,che.dat,geo.dat,5000,1,1);
% 
% ==================================================
%       Variation Partitioning Analysis (VPA):        
% --------------------------------------------------
% 
%     Fraction      R2adj        p   
%     ________    _________    ______
%     'A'           0.03633    0.0554
%     'B'          0.076087    0.0218
%     'C'         0.0006274     0.407
%     'D'           0.10752       NaN
%     'E'           0.16526       NaN
%     'F'          0.075717       NaN
%     'G'           0.12552       NaN
%     'H'           0.41294       NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.345  
%    W = 0.474  
%    Z = 0.367  
% --------------------------------------------------


% Save as PDF:
f_pdf('Doubs_VPA_3')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Distance-based RDA:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpaDB = f_vpaDB(f_dis(bio.dat,'hel'),size(bio.dat,2),top.dat,che.dat,geo.dat,5000,1,1);
% 
% ==================================================
%     Variation Partitioning Analysis (db-VPA):       
% --------------------------------------------------
% 
%     Fraction      R2adj        p   
%     ________    _________    ______
%     'A'           0.03633    0.0574
%     'B'          0.076087    0.0238
%     'C'         0.0006274    0.3872
%     'D'           0.10752       NaN
%     'E'           0.16526       NaN
%     'F'          0.075717       NaN
%     'G'           0.12552       NaN
%     'H'           0.41294       NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.345  
%    W = 0.474  
%    Z = 0.367  
% --------------------------------------------------
