python ./create_batchfolder_lookup_table.py --number_of_batchfolders 10
Start-Process batch_arguments.bat -ArgumentList @(0, 4)
Start-Process batch_arguments.bat -ArgumentList @(1, 4)
Start-Process batch_arguments.bat -ArgumentList @(2, 4)
Start-Process batch_arguments.bat -ArgumentList @(3, 4)


