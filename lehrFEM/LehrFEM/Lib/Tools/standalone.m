function  stand=standalone(main_file, svn_dir, target_dir)
%   Copyright 2006 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%
%standalone(main_file, svn_dir, target_dir) is a script that searches for
%m-files necessary to run the m-file main_file, and copies those files that
%are stored in the repository directory (svn_dir) into a target directory
%target_dir. !!! dat files are taken into account. !!!!!
%standalone('main_QFE','/scratch/users/hheumann/svn/','/scratch/users/hheum
%ann/QFE_auto/')
  

  main_file
  svn_dir
  target_dir
  list=depfun(main_file);
  mkdir(target_dir);
  dep=regexp(list,svn_dir);
  res=find(abs(cellfun('isempty',dep)-1));
  fprintf('Copy the following files to %s', target_dir);
  for i=1:size(res,1);
     list{res(i)}
     copyfile(list{res(i)},target_dir);
  end;
  
return