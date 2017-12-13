main_file=' main_LFE.m'
call=strcat('grep -e [_a-zA-Z0-9]*\.dat -o', main_file) 
[status,result] = system(call)
datfiles=regexp(result,'[_a-zA-Z0-9]*\.dat -o')