$ErrorActionPreference = "Stop"
python -m src.train.train_model_c
python -m src.export.fluent_table
python -m src.export.validate_fluent_table
python -m unittest discover -s tests -v
Write-Host "Pipeline complete. Fluent files are under udf/ and artifacts/fluent_table/."
