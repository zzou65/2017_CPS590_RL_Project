function  stand=create_LehrFEM_Doku(LehrFEMDIR)
%   Create html documentation of LehrFEM
%
%   create_LehrFEM_Doku(LehrFEMDIR) uses m2html (see
%   help m2html) to create documentation
%
%   Example: 
%        create_LehrFEM_Doku('/home/hheumann/svn_local/LehrFEM')
%
%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%   create_LehrFEM_Doku('/home/hheumann/svn_local/LehrFEM')


% change to LehrFEM directory
cd(LehrFEMDIR);
% one more up
cd ..;
% current directory
root=pwd;
% relative LehrFem directory 
LehrFEMDIR = strrep(LehrFEMDIR, [root,'/'], '');
% call m2html
m2html('mfiles',{[LehrFEMDIR,'/Lib'],[LehrFEMDIR,'/Examples']},...
    'htmldir','LehrFEM_DOKU','recursive','on','global','on',...
    'template','frame', 'index','menu');