% This example demonstrates the equivalence of using either traditional or
% distance-based RDA to perform Variation Partitioning Analysis
% 
% by David L. Jones, Apr-2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example follows the analysis presented in Borcard et al. 2011: Section
% 6.3.2 (see p. 183).
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
% Partition the environmental data from the Doubs River into two sets of
% variables: Physiographic (alt, pen, and deb) and Chemical (ph, dur, pho, nit,
% amm, oxy, and dbo).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create three subsets of the env data:
phy.dat = env.dat(:,2:4);  % physiography (upstream-downstream gradient)
phy.txt = env.txt(2:4);

che.dat = env.dat(:,5:11); % chemistry (water quality)
che.txt = env.txt(5:11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Traditional RDA:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hellinger-transform the biotic data:
vpa = f_vpa(f_hellinger(bio.dat),che.dat,phy.dat,[],1000,1,1);
% 
% ==================================================
%       Variation Partitioning Analysis (VPA):        
% --------------------------------------------------
% 
%     Fraction     R2adj       p  
%     ________    _______    _____
%     'A'         0.24135    0.001
%     'B'         0.23304      NaN
%     'C'         0.11205    0.001
%     'D'         0.41356      NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.474  
%    W = 0.345  
% --------------------------------------------------

% Save as PDF:
f_pdf('Doubs_VPA_2')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Distance-based RDA:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpaDB = f_vpaDB(f_dis(bio.dat,'hel'),size(bio.dat,2),che.dat,phy.dat,[],1000,1,1);
% 
% ==================================================
%     Variation Partitioning Analysis (db-VPA):       
% --------------------------------------------------
% 
%     Fraction     R2adj       p  
%     ________    _______    _____
%     'A'         0.24135    0.001
%     'B'         0.23304      NaN
%     'C'         0.11205    0.003
%     'D'         0.41356      NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.474  
%    W = 0.345  
% --------------------------------------------------
