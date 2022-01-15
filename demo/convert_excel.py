import pandas as pd

read_file = pd.io.excel.read_excel (r'frailty/Frailty data from Schultz_Kane_et_al.xlsx')
read_file.to_csv (r'Schultz.csv', index = None, header=True)
