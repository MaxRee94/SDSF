param(
	[Parameter(Mandatory=$true)][string]$config = "default_config.json"
)

python ./create_batchfolder_lookup_table.py --number_of_batchfolders 10
Start-Process batch_arguments.bat -ArgumentList @(0, $config, 4)
Start-Process batch_arguments.bat -ArgumentList @(1, $config, 4)
Start-Process batch_arguments.bat -ArgumentList @(2, $config, 4)
Start-Process batch_arguments.bat -ArgumentList @(3, $config, 4)


