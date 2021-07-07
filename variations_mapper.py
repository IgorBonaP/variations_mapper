#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns

from pathlib import Path

sns.set_palette("Blues")

# File names suffix used by GeneWiz Genome basecalls .csv files.
GENOMIC_VARIATIONS_FILE_NAMES = {
	"SNV": ".Genome_SNV.xls",
	"INDEL": ".Genome_InDel.xls"
}

# Column names used by GeneWiz basecalls .csv files.
GENOMIC_VARIATIONS_COLUMNS = [
	"Type",
	"Chr",
	"Start",
	"END",
	"Ref",
	"Obs",
	"Func",
	"Gene",
	"ExonicFunc",
	"AAChange",
	"hom/het",
	"Qual",
	"Depth",
	"Freq"
]

# Flag to ignore files
flag = "unique_variations"

spreadsheets_folder = r"datasets\basecalls"
output_folder = r"datasets\output"

MAP_RESOLUTION = 1_000_000 # Size of chromosome map bins in megabases.

YAXIS_TICKS = range(0, 10, 2)

CHR_SIZES = {
	# Size is computed taking into account the aproximate size of the chromosome.
	"I": 16_000_000,
	"II": 16_000_000,
	"III": 14_000_000,
	"IV": 18_000_000,
	"V": 21_000_000,
	"X": 18_000_000,
}

CHR_XAXIS_TICKS = {
	chromosome: range(0, size+1, int(size/4))
	for chromosome, size in CHR_SIZES.items()
}

CHR_SIZES_BINS = {chromosome: len(range(0, size, MAP_RESOLUTION))
for chromosome, size in CHR_SIZES.items()}

CHR_SIZES_BINS__ = {chromosome: pd.interval_range(
	start=0,
	end=size,
	freq=MAP_RESOLUTION,) for chromosome, size in CHR_SIZES.items()}
#%%
def load_genome_dataset(strain_name: str, spreadsheets_folder: str)->pd.DataFrame:
	'''
	Load GeneWiz .Genome.csv datasets.
	Indexes the dataset with a multiindex using.
	GENOMIC_VARIATIONS_COLUMNS to uniquely identify each variant:
		"Chr",
		"Type",
		"Start",
		"END",
		"Ref",
		"Obs"
	'''
	snv_ss_path = Path(spreadsheets_folder).joinpath(
		f"{strain_name}{GENOMIC_VARIATIONS_FILE_NAMES['SNV']}")

	indel_ss_path = Path(spreadsheets_folder).joinpath(
		f"{strain_name}{GENOMIC_VARIATIONS_FILE_NAMES['INDEL']}")

	snvs_df = pd.read_csv(snv_ss_path,
		sep="\t",
		usecols=GENOMIC_VARIATIONS_COLUMNS,)
		#index_col=GENOMIC_VARIATIONS_COLUMNS[0:6])

	indel_df = pd.read_csv(indel_ss_path,
		sep="\t",
		usecols=GENOMIC_VARIATIONS_COLUMNS,)
		#index_col=GENOMIC_VARIATIONS_COLUMNS[0:6])

	return snvs_df.append(indel_df, ignore_index=True).set_index(
		keys=[
		"Chr",
		"Type",
		"Start",
		"END",
		"Ref",
		"Obs",]).sort_values(by="Start", axis=0)

def get_EMS_changes(df):
	'''
	Outputs a dataframe containing only the EMS caused SNVS (C/G > T/A)
	'''
	CT_mask = (df.index.get_level_values("Ref") == "C") & (df.index.get_level_values("Obs") == "T")
	GA_mask = (df.index.get_level_values("Ref") == "G") & (df.index.get_level_values("Obs") == "A")
	return df.loc[CT_mask | GA_mask]

def setup_axes_params(ax, title, bg_markers:dict=None)->None:
	ax.set_frame_on(True)
	ax.set_title(title, fontsize=14, weight="bold", family="serif")
	ax.set_yticks(YAXIS_TICKS)
	ax.tick_params(axis="y",
		which= "both",
		length=2,
		direction="in",
		pad=2,
		labelsize="large",
		grid_linestyle="--",
		grid_visible=True,)
	ax.tick_params(axis="x",
		which= "both",
		length=2,
		pad=5,
		rotation=45,
		labeltop=False,
		labelsize="large",
		labelbottom=True,
		grid_linestyle="--",
		grid_visible=True,)
	#ax.ticklabel_format(
	#	axis="x",
	#	style="sci",
	#	scilimits=(0,7))

	#ax.autoscale(axis="y", enable=False)
	ax.autoscale(axis="x", enable=False)
	if bg_markers and bg_markers.get(title, None):
		for name, pos in bg_markers.get(title).items():
			ax.axvline(x=pos,
			label=name,
			color = "black",)
		ax.legend(bbox_to_anchor=(1.05, 1),
		loc='upper left')

def get_map(df:pd.DataFrame,
	strain:str,
	EMS=True,
	hom=False,
	savecountsfolder=False,
	diff_dataset=True,
	bg_markers:dict = None):
	'''
	Returns figure object with the histogram of each chromosome
	in DF.
	'''
	genotype_col_name = strain
	genotype_condition = "het"

	if  not diff_dataset:
		genotype_col_name = GENOMIC_VARIATIONS_COLUMNS[10]
	if hom:
		genotype_condition = "hom"

	if EMS:
		_df = df.pipe(get_EMS_changes).loc[df[genotype_col_name] == genotype_condition]
		variation_type = "EMS"
	else:
		_df = df.loc[df[genotype_col_name] == genotype_condition]
		variation_type = "ALL"

	chr_groups = _df.groupby(level="Chr")
	fig, axes = plt.subplots(
		nrows=1,
		ncols=chr_groups.ngroups,
		figsize=(20, 5),
		sharey=True)

	print(f"Plotting maps of {variation_type} variations.",
	f"Bins width of approx. {MAP_RESOLUTION} Bases.", sep="\n")
	for i, chr_ in enumerate(chr_groups):
		if not isinstance(axes, np.ndarray):
			ax = axes
		else:
			ax = axes[i]

		setup_axes_params(ax, chr_[0], bg_markers=bg_markers)
		count_variations(chr_[1], chr_[0])
		sns.distplot(
		chr_[1].index.get_level_values("Start"),
		axlabel=" ",
		bins=CHR_SIZES_BINS[chr_[0]],
		hist=True,
		vertical=False,
		kde=False,
		norm_hist=False,
		ax=ax)

	return fig

def get_map_fixed_bin(
	df:pd.DataFrame,
	strain:str,
	EMS=True,
	hom=False,
	savecountsfolder=False,
	diff_dataset=True,
	bg_markers:dict = None):
	'''
	Returns figure object with the histogram of each chromosome
	in DF.
	'''
	genotype_col_name = strain
	genotype_condition = "het"

	if  not diff_dataset:
		genotype_col_name = GENOMIC_VARIATIONS_COLUMNS[10]
	if hom:
		genotype_condition = "hom"

	if EMS:
		_df = df.pipe(get_EMS_changes).loc[df[genotype_col_name] == genotype_condition]
		variation_type = "EMS"
	else:
		_df = df.loc[df[genotype_col_name] == genotype_condition]
		variation_type = "ALL"

	chr_groups = _df.groupby(level="Chr")
	fig, axes = plt.subplots(
		nrows=1,
		ncols=chr_groups.ngroups,
		figsize=(40, 5),
		sharey=True,
		tight_layout=True)

	print(f"Plotting maps of {variation_type} variations.")
	for i, chr_ in enumerate(chr_groups):
		if not isinstance(axes, np.ndarray):
			ax = axes
		else:
			ax = axes[i]

		setup_axes_params(ax, chr_[0], bg_markers=bg_markers)
		count_variations(chr_[1], chr_[0])

		get_discrete_histogram(
			chr_[1],
			chr_[0],
			ax=ax)

	return fig

def filter_variations(
	working_strain:pd.DataFrame,
	bg_strain_df:pd.DataFrame,
	operation:str)->pd.DataFrame:
	'''
	Filters working strain basecalls df
	for unique or common mutations with background strain basecalls df.
	'''
	if operation.lower() == "unique":
		mask = ~working_strain.index.isin(bg_strain_df.index)
	elif operation.lower() == "common":
		mask = working_strain.index.isin(bg_strain_df.index)
	else:
		raise NotImplementedError(f"{operation.lower()} is invalid.")
	return working_strain[mask]

def count_variations(df, title, savepath=None):
	'''
	Counts the number of variations in a given basecall df.
	'''
	variations_count = {title: df.index.get_level_values("Start").size}
	print(variations_count)
	return variations_count

def df_dict_to_csv(df_dict:dict, savepath, title_flags=[])->None:
	'''
	Saves all dfs contained within df_dict as csvs into savepath.
	'''
	for title, df in df_dict.items():
		filename = f"{title}.{'.'.join(title_flags)}.csv"
		df.to_csv(Path(savepath).joinpath(filename))

def get_unique_genome_variations(strain: str, bg_strains: list)->dict:
	'''
	Returns df with variations diff between strain
	and all strains in the bg_strains list.
	'''
	strain_df = load_genome_dataset(strain, spreadsheets_folder)

	print("Diffing ", strain)
	print("\tTotal SNVS: ", len(strain_df.index))
	print("\tBackground strains: ", bg_strains)

	for bg_strain in bg_strains:
		count = len(strain_df.index)
		bg_df = load_genome_dataset(bg_strain, spreadsheets_folder)
		common = len(strain_df.pipe(filter_variations, bg_df, "common").index)
		strain_df = strain_df.pipe(filter_variations, bg_df, "unique")
		print(
			f"[{bg_strain}]\t",
			count,
			"-",
			common,
			"=",
			len(strain_df.index))

	print("\tUnique SNVS: ", len(strain_df.index), end="\n")
	return strain_df

def get_last_index_level_value(df:pd.DataFrame, level: int):
	return df.iloc[-1].name[level]

def get_megabase_bin_size(df, megabases=1):
	chr_size = get_last_index_level_value(df, 2)
	bins = range(0,chr_size, np.multiply(megabases, 1000000))
	return bins

def get_discrete_histogram(df, chromosome, ax=None):
	'''
	Returns a histogram binned by chromosome mB with the count
	of variations per mB.
	'''
	variations_pos = df.index.get_level_values("Start").values
	buckets = pd.Series(pd.cut(
		df.index.get_level_values("Start"),
		CHR_SIZES_BINS__[chromosome]))
	cmap = mpl.cm.cool
	frequency = buckets.groupby(buckets).count()
	plot = sns.barplot(
		y=frequency.values,
		x=frequency.index,
		color="slategray",
		ax=ax,)
	xlabels = [np.multiply(n,MAP_RESOLUTION/1_000_000) for n in range(0, len(frequency.index))]
	plot.set_xticklabels(xlabels)
	return plot

def get_shared_mutated_genes(
	working_strain:pd.DataFrame,
	bg_strain_df:pd.DataFrame)->pd.DataFrame:
	shared_df = working_strain.merge(
		bg_strain_df,
		how="inner",
		on=GENOMIC_VARIATIONS_COLUMNS[7],
	)

	return shared_df.loc[(shared_df[GENOMIC_VARIATIONS_COLUMNS[6]] == "exonic")
		& (shared_df[GENOMIC_VARIATIONS_COLUMNS[8]] != "nonsynonymous SNV")
		& (shared_df[GENOMIC_VARIATIONS_COLUMNS[10]] != "hom")]
