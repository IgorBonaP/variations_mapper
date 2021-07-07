#%% [markdown]
#### Imports
#%%
import variations_mapper as vmpr
import pandas as pd
from pathlib import Path
import strains_and_landmarks as SnL
#%% [markdown]
#### Parameters
#%%
strain_of_interest = "QH8150"
df_of_interest = vmpr.load_genome_dataset(strain_of_interest, vmpr.spreadsheets_folder)
savepath = Path(fr"datasets\{strain_of_interest}")

if not savepath.exists():
	savepath.mkdir()
#%%
comparisons_count = []
comparisons_strains = []
for s in SnL.BACKGROUND_STRAINS:
	comparison_strain_df = vmpr.load_genome_dataset(s, vmpr.spreadsheets_folder)
	common = vmpr.filter_variations(df_of_interest, comparison_strain_df, "common")
	comparisons_count.append(len(common.index))
	comparisons_strains.append(s)
comparison_df = pd.DataFrame(comparisons_count, index=comparisons_strains, columns=["COUNT"])
comparison_df.sort_values("COUNT", inplace=True)
comparison_df.plot.bar()

#%% [markdown]
#### Number of common variations between strains...
#%%
comparison_df
#%% [markdown]
#### Filtering unique variations
#%%
unique_variations_strain = vmpr.get_unique_genome_variations(strain_of_interest, SnL.BACKGROUND_STRAINS)
unique_variations_strain.to_csv(savepath.joinpath(f"{strain_of_interest}.unique_variations.csv"))
#%% [markdown]
#### Total unique variations
#%%
len(unique_variations_strain.index)

#%% [markdown]
#### All homozygous nonsynonimous mutations...
#%%
hom_nonsyn_map = vmpr.get_map_fixed_bin(
	unique_variations_strain[unique_variations_strain["ExonicFunc"] != "synonymous SNV"],
	strain_of_interest,
	hom=True,
	EMS=False,
	diff_dataset=False,
	bg_markers=SnL.CHR_LANDMARKS)

#%% [markdown]
#### Homozygous EMS mutations...
#%%
hom_EMS_map = vmpr.get_map_fixed_bin(
	unique_variations_strain,
	strain_of_interest,
	EMS=True,
	hom=True,
	diff_dataset=False,
	bg_markers=SnL.CHR_LANDMARKS)
#%% [markdown]
#### Homozygous nonsynonimous exonic mutations...
#%%
hom_exonic_mutation_map = vmpr.get_map_fixed_bin(
	unique_variations_strain[(unique_variations_strain["Func"] == "exonic") &
	(unique_variations_strain["ExonicFunc"] != "synonymous SNV")],
	strain_of_interest,
	EMS=False,
	hom=True,
	diff_dataset=False,
	bg_markers=SnL.CHR_LANDMARKS)
#%%