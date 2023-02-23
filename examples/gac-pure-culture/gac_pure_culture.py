import pandas as pd
import matplotlib.pyplot as plt

no_gac = pd.read_csv('no_gac.tsv', sep='\t', skiprows=1, usecols=range(1, 11), index_col=0, header=None, names=[
    'step', '# atoms', 'Pressure', 'Mass', '# met', '# eps', 'c_h2', 'c_gco2', 'c_co2', 'c_ch4'])
sheet = pd.read_csv('sheet.tsv', sep='\t', skiprows=1, usecols=range(1, 11), index_col=0, header=None, names=[
    'step', '# atoms', 'Pressure', 'Mass', '# met', '# eps', 'c_h2', 'c_gco2', 'c_co2', 'c_ch4'])
sinusoid = pd.read_csv('sinusoid.tsv', sep='\t', skiprows=1, usecols=range(1, 11), index_col=0, header=None, names=[
    'step', '# atoms', 'Pressure', 'Mass', '# met', '# eps', 'c_h2', 'c_gco2', 'c_co2', 'c_ch4'])
i = max(len(no_gac), len(sheet), len(sinusoid))

if len(no_gac) < i:
    no_gac = pd.concat([no_gac, pd.DataFrame([no_gac.iloc[-1].values.tolist()] * (i - len(no_gac)), columns=no_gac.columns, index=range(len(no_gac), i))])
if len(sheet) < i:
    sheet = pd.concat([sheet, pd.DataFrame([sheet.iloc[-1].values.tolist()] * (i - len(sheet)), columns=sheet.columns, index=range(len(sheet), i))])
if len(sinusoid) < i:
    sinusoid = pd.concat([sinusoid, pd.DataFrame([sinusoid.iloc[-1].values.tolist()] * (i - len(sinusoid)), columns=sinusoid.columns, index=range(len(sinusoid), i))])


# Methane total
x = sheet.index
no_gac_ch4 = no_gac['c_ch4']
sheet_ch4 = sheet['c_ch4']
sinusoid_ch4 = sinusoid['c_ch4']
plt.clf()
plt.plot(x, no_gac_ch4, label='No GAC')
plt.plot(x, sheet_ch4, label='GAC sheet')
plt.plot(x, sinusoid_ch4, label='GAC sinusoid')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Methane Concentration (g/m3)')
plt.title('Methane Concentration over Time')
plt.savefig('methane_by_gac_shape.png', dpi=300, bbox_inches='tight')
plt.clf()

# Methane per cell
no_gac_ch4 = no_gac['c_ch4'] / no_gac['# met']
sheet_ch4 = sheet['c_ch4'] / sheet['# met']
sinusoid_ch4 = sinusoid['c_ch4'] / sinusoid['# met']
plt.plot(x, no_gac_ch4, label='No GAC')
plt.plot(x, sheet_ch4, label='GAC sheet')
plt.plot(x, sinusoid_ch4, label='GAC sinusoid')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Methane Concentration per Cell (g/m3 * #met-1)')
plt.title('Methane Concentration over Time per Number of cells')
plt.savefig('methane_by_number_cells.png', dpi=300, bbox_inches='tight')
