python ./create_batchfolder_lookup_table.py --number_of_batchfolders 10
Start-Process batch_arguments.bat -ArgumentList @(0, 8) -RedirectStandardOutput "log_proc0.txt" -RedirectStandardError "logError_proc0.txt"
Start-Process batch_arguments.bat -ArgumentList @(1, 8) -RedirectStandardOutput "log_proc1.txt" -RedirectStandardError "logError_proc1.txt"
Start-Process batch_arguments.bat -ArgumentList @(2, 8) -RedirectStandardOutput "log_proc2.txt" -RedirectStandardError "logError_proc2.txt"
Start-Process batch_arguments.bat -ArgumentList @(3, 8) -RedirectStandardOutput "log_proc3.txt" -RedirectStandardError "logError_proc3.txt"
Start-Process batch_arguments.bat -ArgumentList @(4, 8) -RedirectStandardOutput "log_proc4.txt" -RedirectStandardError "logError_proc4.txt"
Start-Process batch_arguments.bat -ArgumentList @(5, 8) -RedirectStandardOutput "log_proc5.txt" -RedirectStandardError "logError_proc5.txt"
Start-Process batch_arguments.bat -ArgumentList @(6, 8) -RedirectStandardOutput "log_proc6.txt" -RedirectStandardError "logError_proc6.txt"
Start-Process batch_arguments.bat -ArgumentList @(7, 8) -RedirectStandardOutput "log_proc7.txt" -RedirectStandardError "logError_proc7.txt"


