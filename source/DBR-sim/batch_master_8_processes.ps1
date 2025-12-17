.\env\Scripts\Activate.ps1
python ./create_batchfolder_lookup_table.py --number_of_batchfolders 10
Start-Process batch_arguments.bat -ArgumentList @(0, 8)
Start-Process batch_arguments.bat -ArgumentList @(1, 8)
Start-Process batch_arguments.bat -ArgumentList @(2, 8)
Start-Process batch_arguments.bat -ArgumentList @(3, 8)
Start-Process batch_arguments.bat -ArgumentList @(4, 8)
Start-Process batch_arguments.bat -ArgumentList @(5, 8)
Start-Process batch_arguments.bat -ArgumentList @(6, 8)
Start-Process batch_arguments.bat -ArgumentList @(7, 8)


