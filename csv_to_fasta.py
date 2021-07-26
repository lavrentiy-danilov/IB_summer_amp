import pandas as pd


df = pd.read_csv('general_amps-1.csv', header=None, skiprows=[0])
df = df.iloc[:, 0:2]
df[0] = df[0].apply(lambda x: '>' + x)
df[1] = df[1].apply(lambda x: x.upper().replace('"', '').replace('\n', ''))
m = df[1].str.match(r'[A-Z]+\b')
df = df[m]
print(m.head())
print(df.head())

df.to_csv('database.fasta', sep='\n', index=False, header=False)
