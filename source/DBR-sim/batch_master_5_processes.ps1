python ./create_batchfolder_lookup_table.py --number_of_batchfolders 10
Start-Process batch_arguments.bat -ArgumentList @(0, 5)
Start-Process batch_arguments.bat -ArgumentList @(1, 5)
Start-Process batch_arguments.bat -ArgumentList @(2, 5)
Start-Process batch_arguments.bat -ArgumentList @(3, 5)
Start-Process batch_arguments.bat -ArgumentList @(4, 5)


