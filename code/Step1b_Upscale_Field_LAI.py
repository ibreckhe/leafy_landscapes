# Import necessary libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gspread_pandas import Spread, Client
from oauth2client.service_account import ServiceAccountCredentials

# Set working directory
wd = "/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents\ -\ Research\ -\ Spatial\ Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/"

# Read leaf area data and plot data from CSV and Google Sheets, respectively
area_df_sum = pd.read_csv(wd + "data/plot_leaf_areas_8_8_2024.csv")

# Authenticate with Google Sheets
scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('path_to_your_service_account_file.json', scope)
client = Client(scope=scope, creds=creds)
spread = Spread('1rSPt5g_-af3EsmYj_nXwJO9Ya9uL4duAkrBlrsVNMSo')
plot_data = spread.sheet_to_df(sheet="PL_data", index=None, start_row=4, end_row=317)

# Format plot data for joining
plot_data['PlotID'] = plot_data['PlotID'].str.upper()
plot_data['Num_leaves_scanned'] = pd.to_numeric(plot_data['Num_leaves_scanned'])

# Join leaf area data with plot data
plot_data_area = pd.merge(plot_data[['PlotID', 'Species_or_FT', 'Collected_date', 'Num_leaves_scanned', 'Scanned_wet_mass', 'Scanned_dry_mass', 'Scanned_leaf_water_pct']],
                          area_df_sum.rename(columns={'plot': 'PlotID', 'sample': 'Species_or_FT'}),
                          on=['PlotID', 'Species_or_FT'],
                          how='left')

# Calculate LMA for each plot
plot_data_area['LMA'] = plot_data_area['Scanned_dry_mass'] / (plot_data_area['total_leaf_area_sqcm'] * 0.0001)

# Plot LMA for each functional type
sns.set_style("whitegrid")
p1 = sns.stripplot(x='Species_or_FT', y='LMA', hue='Species_or_FT', data=plot_data_area[plot_data_area['Species_or_FT'].isin(['forb', 'grass', 'shrub'])],
                   dodge=True, alpha=0.1)
p1.set_title("Functional Group LMA")
p1.set_ylabel("CWM Leaf Mass Per Area [g/m^2]")
p1.set_ylim(0, 420)

# Plot LMA vs. leaf water content
p2 = sns.lmplot(x='Scanned_leaf_water_pct', y='LMA', hue='Species_or_FT', data=plot_data_area[plot_data_area['Species_or_FT'].isin(['forb', 'grass', 'shrub'])],
                 scatter_kws={'alpha':0.1})
p2.set(xlabel="CWM Leaf Water Content [%]", ylabel="CWM Leaf Mass Per Area [g/m^2]", ylim=(0, 420))
plt.suptitle("LMA vs. Leaf Water Content")

# Save plots as PDF
plt.savefig(wd + "../LMA_plots.pdf")

# Write joined data to disk
plot_data_area.to_csv(wd + "data/finalized_plot_data/leafy_landscapes_plot_data_v1.csv", index=False)