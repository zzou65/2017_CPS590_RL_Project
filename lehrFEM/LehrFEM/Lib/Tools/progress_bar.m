function [] = progress_bar(i)
% PROGRESS_BAR Displays a progress bar in the command window.
%
%   PROGRESS_BAR(PER) displays the progress percentage PER as a progress
%   bar inside the command window.
%
%   Example:
%
%   progress_bar(0);
%   for i = 1:100
%     progress_bar(i);
%     pause(.1);
%   end

%   Copyright 2006-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Clear previous progress bar from command window

  if(i >= 1)  
    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' ...
             '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b']);
  end
  
  % Display current progress
  
  if(0 <= i && i < 5)
    fprintf('[                                        ]  :  %3d %%', i); 
  elseif(5 <= i && i < 10)
    fprintf('[##                                      ]  :  %3d %%', i);   
  elseif(10 <= i && i < 15)
    fprintf('[####                                    ]  :  %3d %%', i);    
  elseif(15 <= i && i < 20)
    fprintf('[######                                  ]  :  %3d %%', i);        
  elseif(20 <= i && i < 25)
    fprintf('[########                                ]  :  %3d %%', i);    
  elseif(25 <= i && i < 30)
    fprintf('[##########                              ]  :  %3d %%', i);  
  elseif(30 <= i && i < 35)
    fprintf('[############                            ]  :  %3d %%', i);  
  elseif(35 <= i && i < 40)
    fprintf('[##############                          ]  :  %3d %%', i);  
  elseif(40 <= i && i < 45)
    fprintf('[################                        ]  :  %3d %%', i);  
  elseif(45 <= i && i < 50)
    fprintf('[##################                      ]  :  %3d %%', i);  
  elseif(50 <= i && i < 55)
    fprintf('[####################                    ]  :  %3d %%', i);  
  elseif(55 <= i && i < 60)
    fprintf('[######################                  ]  :  %3d %%', i);  
  elseif(60 <= i && i < 65)
    fprintf('[########################                ]  :  %3d %%', i);  
  elseif(65 <= i && i < 70)
    fprintf('[##########################              ]  :  %3d %%', i);  
  elseif(70 <= i && i < 75)
    fprintf('[############################            ]  :  %3d %%', i);
  elseif(75 <= i && i < 80)
    fprintf('[##############################          ]  :  %3d %%', i);  
  elseif(80 <= i && i < 85)
    fprintf('[################################        ]  :  %3d %%', i);  
  elseif(85 <= i && i < 90)
    fprintf('[##################################      ]  :  %3d %%', i);  
  elseif(90 <= i && i < 95)
    fprintf('[####################################    ]  :  %3d %%', i);  
  elseif(95 <= i && i < 100)
    fprintf('[######################################  ]  :  %3d %%', i);  
  else    
    fprintf('[########################################]  :  100 %%');
  end

  % Delete progress bar from command window
  
  if(i == 100)
    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' ...
             '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b']);
  end
    
return