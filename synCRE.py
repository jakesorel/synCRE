import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys,os
import json
import gzip, shutil
from sklearn.cluster import KMeans
from io import StringIO
from pybedtools import BedTool
import multiprocessing
from joblib import Parallel, delayed
import seaborn as sb

class Lookup:
    """

    """
    def __init__(self,gene_info="reference/gene_info",
                 TF_list="reference/TF_list.csv",
                 RNA_seq_file="reference/RNA_seq.txt",
                 motif_annotations="reference/motif_annotations.xlsx",
                 manual_addition_file="reference/manual_additions.csv"):

        self.gene_info = gene_info
        self.TF_list = TF_list
        self.RNA_seq_file = RNA_seq_file
        self.motif_annotations = motif_annotations
        self.manual_addition_file = manual_addition_file

        make_directory("reference/lookup_table")

        self.all_gene_list = pd.read_csv(self.gene_info, sep="\t")
        self.true_names = capitalize_list(list(self.all_gene_list["Symbol"]))
        self.make_alias_dict()

        self.make_lookup()
        self.save_lookup()
        self.save_not_referenced()

    def make_alias_dict(self):
        """

        :return:
        """
        self.find_all_names()
        self.alias_dict_generator()
        self.save_aliases()

    def find_all_names(self):
        """

        :return:
        """
        self.all_names = []
        for i, synonym in enumerate(self.all_gene_list["Synonyms"].values):
            if synonym == "-":
                self.all_names.append([self.true_names[i]])
            else:
                syn_vals = capitalize_list(synonym.split("|"))
                self.all_names.append([self.true_names[i]] + syn_vals)

    def alias_dict_generator(self):
        """

        :return:
        """
        self.aliases = {}
        for i, names in enumerate(self.all_names):
            for name in names:
                self.aliases[name] = self.true_names[i]


    def save_aliases(self):
        """

        :return:
        """
        save_dict(self.aliases, "reference/lookup_table/aliases")

    def load_TF_names(self):
        """

        :return:
        """
        self.TF_names = list(pd.read_csv(self.TF_list,header=None)[0].values)
        self.TF_names = self.standardize_names(self.TF_names)

    def load_RNA_names(self):
        """

        :return:
        """
        self.RNA_names = list(pd.read_csv(open(self.RNA_seq_file).read(), index_col=0).index)
        self.RNA_names = self.standardize_names(self.RNA_names)

    def standardize_names(self,names):
        """

        :param names:
        :return:
        """
        return standardize_names(names, self.aliases, self.true_names)

    def load_motif_annotations(self):
        """

        :return:
        """
        self.motifs = pd.read_excel(self.motif_annotations, 1, engine='openpyxl')

    def make_lookup(self):
        """

        :return:
        """
        self.load_motif_annotations()
        self.load_TF_names()
        self.load_RNA_names()
        self.cluster_ids = self.motifs["Cluster_ID"].values
        self.motif_names = self.motifs["Motif"]

        self.split_doubles()
        self.reformat_names()
        self.fix_Znf()
        self.motif_names = self.standardize_names(self.motif_names)
        self.fix_manual_additions()
        self.de_duplicate()


    def split_doubles(self):
        """

        :return:
        """
        motif_names_old, cluster_ids_old = self.motif_names.copy(), self.cluster_ids.copy()
        motif_names, cluster_ids = [], []
        for i, name in enumerate(motif_names_old):
            if len(name.split("+")) > 1:
                names = name.split("+")
                for name in names:
                    motif_names.append(name)
                    cluster_ids.append(cluster_ids_old[i])
            else:
                motif_names.append(name)
                cluster_ids.append(cluster_ids_old[i])
        self.motif_names,self.cluster_ids = motif_names,cluster_ids

    def fix_Znf(self):
        """

        :return:
        """
        for i, motif in enumerate(self.motif_names):
            if ("Zn" in motif) & ("Znf" not in motif):
                self.motif_names[i] = "Znf" + motif.split("Zn")[1]

    def fix_manual_additions(self):
        """

        :return:
        """
        manual_df = pd.read_csv(self.manual_addition_file, index_col=0)
        manual_dict = dict(zip(manual_df.Name, manual_df.New_name))
        orig_names = manual_df.Name.values
        for i, motif in enumerate(self.motif_names):
            if motif in orig_names:
                self.motif_names[i] = manual_dict[motif]

    def de_duplicate(self):
        """

        :return:
        """
        motif_names_old, cluster_ids_old = np.array(self.motif_names), np.array(self.cluster_ids)
        motif_names, cluster_ids = [], []
        for cluster_id in np.unique(cluster_ids_old):
            id_mask = cluster_id == cluster_ids_old
            unique_names = np.unique(motif_names_old[id_mask])
            for name in unique_names:
                motif_names.append(name)
                cluster_ids.append(cluster_id)
        self.motif_names,self.cluster_ids = motif_names,cluster_ids

    def reformat_names(self):
        """

        :return:
        """
        for i, name in enumerate(self.motif_names):
            self.motif_names[i] = name.split("_")[0]
            self.motif_names[i] = self.motif_names[i].split(".")[0].capitalize()

    def save_lookup(self):
        """

        :return:
        """
        save_dict(dict(zip(self.motif_names, [int(cid) for cid in self.cluster_ids])), "reference/lookup_table/lookup_table")

    def save_not_referenced(self):
        """

        :return:
        """
        make_directory("reference/lookup_table/not_referenced")
        save_csv(sorted(list(set(self.motif_names).difference(set(self.RNA_names)))),"reference/lookup_table/not_referenced/motifs_not_in_RNA.csv")
        save_csv(sorted(list(set(self.motif_names).difference(set(self.TF_names)))),"reference/lookup_table/not_referenced/motifs_not_in_TF_list.csv")
        save_csv(sorted(list(set(self.TF_names).difference(set(self.motif_names)))),"reference/lookup_table/not_referenced/TFs_not_in_motif_list.csv")


class Expression:
    """

    """
    def __init__(self,RNA_seq_file="reference/RNA_seq.txt",n_clusters=8,
                 candidate_genes=["Nkx2-2", "Olig2", "Pax6", "Irx3", "Gli3", "Sox2", "Cdx1", "Cdx2"]):
        self.RNA_seq_file = RNA_seq_file
        self.n_clusters = n_clusters
        self.candidate_genes = candidate_genes

        ##Make directory structure
        make_directory("results")
        make_directory("results/expression")
        make_directory("results/expression/plots")
        make_directory("results/expression/clusters")
        make_directory("results/expression/all_genes")
        make_directory("results/expression/archetypes")

        ##Load lookup table
        self.aliases = open_dict("reference/lookup_table/aliases")
        self.true_names = self.aliases.values()
        self.lookup = open_dict("reference/lookup_table/lookup_table")
        self.gene_names = list(self.lookup.keys())

        ##Load RNA_seq_file
        self.RNA_df = {}
        self.RNA_df["all"] = pd.read_csv(open(self.RNA_seq_file).read(), index_col=0)
        self.RNA_df["all"].index = self.standardize_names(self.RNA_df["all"].index)

    def run_all(self):
        """

        :return:
        """
        self.build_z()
        self.calculate_k_means()
        self.save_outputs()
        self.plot_candidate_genes()
        self.plot_clusters()

    def build_z(self,percentile=25,fc_thresh=0.25):
        """

        :param percentile:
        :param fc_thresh:
        :return:
        """
        self.filter_by_lookup()
        self.make_RNA_matrix()
        self.make_mean_matrix()
        self.filter_by_max_expression(percentile)
        self.filter_by_fold_change(fc_thresh)
        self.calculate_z_score()



    def filter_by_lookup(self):
        """
        filter by genes represented in the lookup table
        :return:
        """
        genes = list(self.RNA_df["all"].index.values)
        TF_mask = [gene in self.gene_names for gene in genes]
        self.RNA_df["lookup"] = self.RNA_df["all"].iloc[np.nonzero(TF_mask)]
        self.TF_genes = {}
        self.TF_genes["all"] = np.array(self.RNA_df["lookup"].index.values)

    def make_RNA_matrix(self):
        """
        Makes a matrix from the RNAseq dataframe
        :return:
        """
        conditions = list(self.RNA_df["lookup"].columns)
        condition_split = np.array([condition.split(".") for condition in conditions])
        condition_split[:, 0] = np.array([condition.split("D")[1] for condition in condition_split[:, 0]])
        condition_split[:, 3] = np.array([condition.split("R")[1] for condition in condition_split[:, 3]])
        condition_split[condition_split[:, 2] == "NMP", 2], condition_split[condition_split[:, 2] == "3", 2], \
        condition_split[condition_split[:, 2] == "M", 2], condition_split[condition_split[:, 2] == "2", 2], \
        condition_split[condition_split[:, 2] == "1", 2] = 0, 10, 20, 30, 40
        condition_split = condition_split.astype(np.float)
        condition_split = np.column_stack(
            (condition_split[:, 2], condition_split[:, 0], condition_split[:, 1], condition_split[:, 3]))

        pos_mat, day_mat, conc_mat, rep_mat = np.meshgrid(np.array((0, 10, 20, 30, 40)), np.array((3, 4, 5, 6)),
                                                          np.array((0, 10, 100, 500)), np.array((1, 2, 3)),
                                                          indexing="ij")
        self.expr_mat = np.ones((pos_mat.shape[0], pos_mat.shape[1], pos_mat.shape[2], pos_mat.shape[3], self.RNA_df["lookup"].shape[0]),
                           dtype=np.float64) * np.nan
        pos_id, day_id, conc_id, rep_id = np.zeros(condition_split.shape[0], dtype=np.int64), np.zeros(
            condition_split.shape[0], dtype=np.int64), np.zeros(condition_split.shape[0], dtype=np.int64), np.zeros(
            condition_split.shape[0], dtype=np.int64)
        for i in range(condition_split.shape[0]):
            a, b, c, d = np.nonzero((condition_split[i, 0] == pos_mat) * (condition_split[i, 1] == day_mat) * (
                        conc_mat == condition_split[i, 2]) * (rep_mat == condition_split[i, 3]))
            pos_id[i] = a[0]
            day_id[i] = b[0]
            conc_id[i] = c[0]
            rep_id[i] = d[0]
            self.expr_mat[a[0], b[0], c[0], d[0]] = self.RNA_df["lookup"][conditions[i]]

    def make_mean_matrix(self):
        """

        :return:
        """
        self.mean_expr_mat = {}
        self.mean_expr_mat["all"] = np.nanmean(self.expr_mat, axis=(2, 3))

    def filter_by_max_expression(self,percentile=25):
        """

        :param percentile:
        :return:
        """
        thresh = np.nanpercentile(self.mean_expr_mat["all"][~np.isnan(self.mean_expr_mat["all"])], percentile)
        print("Min expression threshold = ", thresh)
        expressed_mask = (self.mean_expr_mat["all"] > percentile).any(axis=(0, 1))
        self.mean_expr_mat["expr_filtered"] = self.mean_expr_mat["all"][:, :, expressed_mask]
        self.TF_genes["expr_filtered"] = self.TF_genes["all"][expressed_mask]

    def filter_by_fold_change(self,fc_thresh=0.25):
        """

        :param fc_thresh:
        :return:
        """
        fold_change = np.nanmin(self.mean_expr_mat["expr_filtered"], axis=(0, 1)) / np.nanmax(self.mean_expr_mat["expr_filtered"], axis=(0, 1))
        changed_mask = fold_change < fc_thresh
        self.mean_expr_mat["expr+fc_filtered"] = self.mean_expr_mat["expr_filtered"][:, :, changed_mask]
        self.TF_genes["expr+fc_filtered"] = self.TF_genes["expr_filtered"][changed_mask]

    def calculate_z_score(self):
        """

        :return:
        """
        self.mean_expr_z = {}
        self.mean_expr_z["expr+fc_filtered"] = np.dstack([self.z_score(self.mean_expr_mat["expr+fc_filtered"][:, :, i]) for i in range(self.mean_expr_mat["expr+fc_filtered"].shape[-1])])

    def calculate_k_means(self):
        """

        :return:
        """
        empty_mask = ~np.isnan(self.mean_expr_z["expr+fc_filtered"].reshape(-1, self.mean_expr_z["expr+fc_filtered"].shape[-1])).all(axis=1)
        flat_expr_mat_z = self.mean_expr_z["expr+fc_filtered"].reshape(-1, self.mean_expr_z["expr+fc_filtered"].shape[-1])[empty_mask].T
        self.kmeans = KMeans(n_clusters=self.n_clusters, random_state=0).fit(flat_expr_mat_z)
        self.gene_states = self.kmeans.labels_

    def save_outputs(self):
        """

        :return:
        """
        self.save_all_genes()
        self.save_genes_by_cluster()
        self.save_archetypes_by_cluster()

    def save_all_genes(self):
        """

        :return:
        """
        save_csv(list(self.TF_genes["expr+fc_filtered"]),
                 "results/expression/all_genes/filtered_genes.txt", header=False, index=False)

    def save_genes_by_cluster(self):
        """

        :return:
        """
        for i in np.unique(self.gene_states):
            save_csv(list(self.TF_genes["expr+fc_filtered"][self.gene_states == i]),
                     "results/expression/clusters/cluster %d.txt" % i, header=False, index=False)

    def save_archetypes_by_cluster(self):
        """

        :return:
        """
        for i in np.unique(self.gene_states):
            gene_list = list(self.TF_genes["expr+fc_filtered"][self.gene_states == i])
            cid_list = []
            for gene in gene_list:
                try:
                    cid_list.append(self.lookup[gene])
                except:
                    a = 1
            save_csv(np.array(np.unique(cid_list)), "results/expression/archetypes/archetypes_for_cluster_%d.txt" % i,
                     header=False, index=False)

    def plot_candidate_genes(self,cmap=plt.cm.viridis, vmin=-2, vmax=2):
        """

        :param cmap:
        :param vmin:
        :param vmax:
        :return:
        """
        fig, ax = plt.subplots(4, 2, sharey=True, sharex=True)
        ax = ax.ravel()
        for i, gene in enumerate(self.candidate_genes):
            zmap = self.mean_expr_z["expr+fc_filtered"][:,:,np.nonzero(self.TF_genes["expr+fc_filtered"]==gene)[0][0]]
            self.plot_expression_profile(zmap, gene, ax[i], cmap=cmap, vmin=vmin, vmax=vmax)
        fig.subplots_adjust(hspace=0.5, wspace=0)
        fig.savefig("results/expression/plots/candidate_genes.pdf")

    def plot_clusters(self,cmap=plt.cm.viridis, vmin=-2, vmax=2):
        """

        :param cmap:
        :param vmin:
        :param vmax:
        :return:
        """
        fig, ax = plt.subplots(4, 2, sharex=True, sharey=True)
        ax = ax.ravel()
        for j, i in enumerate(np.unique(self.gene_states)):
            av_expr = np.nanmean(self.mean_expr_z["expr+fc_filtered"][:, :, self.gene_states == i], axis=-1)
            self.plot_expression_profile(av_expr, "Cluster %d" % i, ax[j],cmap=cmap,vmin=vmin,vmax=vmax)
        fig.subplots_adjust(hspace=0.5, wspace=0)
        fig.savefig("results/expression/plots/expression_clusters.pdf")

    def standardize_names(self,names):
        """

        :param names:
        :return:
        """
        names = standardize_names(names, self.aliases, self.true_names)
        return names

    def z_score(self,expr):
        """

        :param expr:
        :return:
        """
        return (expr - np.nanmean(expr)) / np.nanstd(expr)

    def plot_expression_profile(self,zmap, title, ax, cmap=plt.cm.viridis, vmin=-2, vmax=2):
        """

        :param zmap:
        :param title:
        :param ax:
        :param cmap:
        :param vmin:
        :param vmax:
        :return:
        """
        ax.set_title(title)
        ax.imshow(np.flip(zmap, axis=0), vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_yticks(np.arange(5))
        ax.set_yticklabels(np.flip(np.array(["NMP", "p3", "pMN", "p2", "p1"])))
        ax.set_xticks(np.arange(4))
        ax.set_xticklabels(["D3", "D4", "D5", "D6"], rotation=90)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmax=vmax, vmin=vmin))
        sm._A = []
        cl = plt.colorbar(sm, ax=ax, pad=0.05, fraction=0.085, aspect=10, orientation="vertical")
        cl.set_label("z-score")

    def kmeans(self,expr_mat_z):
        """

        :param expr_mat_z:
        :return:
        """
        empty_mask = ~np.isnan(expr_mat_z.reshape(-1, expr_mat_z.shape[-1])).all(axis=1)
        flat_expr_mat_z = expr_mat_z.reshape(-1, expr_mat_z.shape[-1])[empty_mask].T
        kmeans = KMeans(n_clusters=8, random_state=0).fit(flat_expr_mat_z)
        gene_states = kmeans.labels_
        return kmeans, gene_states

class Motif_Finder:
    """

    """
    def __init__(self,
                 genome_dir="reference/genome_dir.txt",
                 motif_annotations="reference/motif_annotations.xlsx",
                 chip_truth="reference/chip_truth.txt"):
        self.genome_dir = open(genome_dir).read()
        self.motif_table = pd.read_excel(motif_annotations, 1, engine='openpyxl')
        self.lookup = open_dict("reference/lookup_table/lookup_table")
        self.chip_truth = pd.read_csv(chip_truth,sep="\t",header=None)
        self.hit_thresh = []

    def make_pmf(self):
        """
        This script parses the set of ".meme" files in the meme_files directory.

        It takes each pmf matrix and saves it to a separate file in the "pmf" directory. Filenames are assigned by the motif-name.

        :return:
        """
        meme_files = os.listdir("reference/motifs/meme_files")
        for memefile in meme_files:
            pmf_on = False
            for i, line in enumerate(open("reference/motifs/meme_files/%s" % memefile).readlines()):
                if "MOTIF" in line:
                    motif_name = line.split("MOTIF ")[1].split("\n")[0]
                    motifs = []
                if (("URL" in line) or (line == "\n")) and (pmf_on == True):
                    pmf_on = False
                    df = pd.DataFrame(motifs).transpose()
                    df.to_csv("reference/motifs/pmf/%s.pmf" % motif_name, sep="\t", header=False, index=False)
                if pmf_on is True:
                    motif = pd.read_csv(StringIO(line), delim_whitespace=True, header=None).values.ravel()
                    motifs.append(motif)
                if "letter-probability matrix" in line:
                    pmf_on = True

    def sample_eCRE_sequence(self):
        """

        :return:
        """
        make_directory("results/fasta")
        make_directory("results/fasta/by_eCRE")
        bed_files = os.listdir("reference/eCRE_locs")
        delete_directory("results/fasta/scrap")
        for bed_file in bed_files:
            make_directory("results/fasta/scrap")
            eCRE_name = bed_file.split(".bed")[0]
            eCRE = BedTool("reference/eCRE_locs/%s" % bed_file)
            chrom = eCRE.to_dataframe()["chrom"].values[0]
            chrfa = "%s/%s.fa.gz" % (self.genome_dir,chrom)
            unzip(chrfa, "results/fasta/scrap/out.fa")
            eCRE.sequence("results/fasta/scrap/out.fa", fo="results/fasta/by_eCRE/%s.fa" % eCRE_name)
            delete_directory("results/fasta/scrap")
            print("fasta extraction for ",eCRE_name, " complete")

    def find_motifs(self,minthresh=4):
        """

        Uses MOODS to find all motif matches above a threshold log-odds score (against pre-defined background rates)

        :param p_vals:
        :return:
        """
        make_directory("results/motifs")
        make_directory("results/motifs/raw")

        runline = """
eval "$(conda shell.bash hook)" \n
source activate moods \n
cd pkgs/moods/scripts \n
pwd \n
python2 moods-dna.py  \
--sep ";" -s ../../../results/fasta/by_eCRE/%s.fa --threshold %.6f \
--lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 \
-m ../../../reference/motifs/pmf/* -o ../../../results/motifs/raw/%s.csv
        """

        eCRE_names = [name.split(".bed")[0] for name in os.listdir("reference/eCRE_locs")]
        for name in eCRE_names:
            os.system(runline % (name, minthresh, name))
            print("Motifs identified for the %s eCRE above thresh=%.6f" % (name, minthresh))

    def plot_motif_distributions(self):

        fig2, ax2 = plt.subplots(figsize=(4, 3))
        raw_files = os.listdir("results/motifs/raw")
        for file in raw_files:
            if ".csv" in file:
                df = pd.read_csv("results/motifs/raw/%s" % file, sep=";", header=None)
                make_directory("results/motifs/hit_score")
                fig, ax = plt.subplots(figsize=(4, 3))
                sb.distplot(df[4], bins=200, ax=ax)
                sb.distplot(df[4], bins=200, ax=ax2, label=file.split(".csv")[0])
                ylim = ax.get_ylim()
                if type(self.hit_thresh) is not list:
                    ax.plot((self.hit_thresh,self.hit_thresh),(0,ylim[1]*2),linestyle=":",color="grey",zorder=10)
                ax.set(xlabel="Hit Score", ylabel="Density",ylim=ylim)

                fig.tight_layout()
                fig.savefig("results/motifs/hit_score/%s.pdf" % file.split(".csv")[0])
        ylim2 = ax2.get_ylim()
        if type(self.hit_thresh) is not list:
            ax2.plot((self.hit_thresh, self.hit_thresh), (0, ylim2[1] * 2), linestyle=":", color="grey", zorder=10)
        ax2.set(xlabel="Hit Score", ylabel="Density", ylim=ylim2)
        ax2.legend()
        fig2.tight_layout()
        fig2.savefig("results/motifs/hit_score/merge.pdf")
    #
    # def make_truth_matrix(self):
    #     eCRE_beds = []
    #     for bed in os.listdir("reference/eCRE_locs"):
    #         if ".DS" not in bed:
    #             eCRE_beds.append(bed)
    #     eCRE_names = [bed.split(".bed")[0] for bed in eCRE_beds]
    #     dicts = []
    #     for bed in eCRE_beds:
    #         names, hits = [],[]
    #         for i in range(self.chip_truth.shape[0]):
    #             chip_name, publication, file = self.chip_truth.iloc[i]
    #             hit = BedTool("reference/eCRE_locs/%s"%bed).intersect(file).to_dataframe()
    #             print(hit)
    #             names.append(chip_name)
    #             hits.append(hit)
    #         dicts.append(dict(zip(names,hits)))
    #     truth_mat = dict(zip(eCRE_names,dicts))
    #     # print(truth_mat)

    def get_threshold(self,required_dicts=False,percentile=80):
        if required_dicts is not False:
            lookup_names = dict(zip(self.motif_table["Motif"], self.motif_table["Cluster_ID"]))
            motifcsvs = os.listdir("results/motifs/raw")
            make_directory("results/motifs/bed")
            bed_dfs = {}
            for motifcsv in motifcsvs:
                if ".csv" in motifcsv:
                    df = pd.read_csv("results/motifs/raw/%s" % motifcsv, sep=";", header=None)
                    seq_len = [len(string) for string in df[5]]
                    locs = df[0]
                    names = df[1].values
                    names = [name.split(".pmf")[0].split(" ")[0] for name in names]
                    motif_ids = []
                    for name in names:
                        try:
                            motif_ids.append(lookup_names[name])
                        except KeyError:
                            motif_ids.append(-1)
                    chrom, start, end = [], [], []
                    for loc in locs:
                        chrm, pos = loc.split(":")
                        st, en = pos.split("-")
                        st, en = int(st), int(en)
                        chrom.append(chrm)
                        start.append(st)
                        end.append(en)
                    start, end = np.array(start), np.array(end)
                    start += df[2].values
                    end = start + seq_len
                    bed_df = pd.DataFrame([chrom, start, end, motif_ids, df[4]]).transpose()
                    bed_dfs[motifcsv.split(".csv")[0]] = bed_df

            thresh_dicts = []
            for i, (eCRE, bed_df) in enumerate(bed_dfs.items()):
                for j, (TF, num) in enumerate(required_dicts[eCRE].items()):
                    a = 1
                    cid = self.lookup[TF]
                    cid_df = bed_df.loc[bed_df[3] == cid]
                    if cid_df.size!=0:
                        thresh_dicts.append(
                            {"eCRE": eCRE, "TF": TF, "thresh": cid_df[4].values[(-1 * cid_df[4]).argsort()][num - 1]})
                    if cid_df.size==0:
                        print("No match found for TF %s in eCRE %s, ignoring..."%(TF,eCRE))
            self.hit_thresh = pd.DataFrame(thresh_dicts)["thresh"].min()
        else:
            raw_files = os.listdir("results/motifs/raw")
            hscores = np.array(())
            for file in raw_files:
                if ".csv" in file:
                    df = pd.read_csv("results/motifs/raw/%s" % file, sep=";", header=None)
                    hscores = np.concatenate((hscores, df[4]))
            self.hit_thresh = np.percentile(hscores, percentile)
        print("thresh = ",self.hit_thresh)

    def filter_motifs_by_hit(self,required_dicts=False,percentile=80):
        self.get_threshold(required_dicts=required_dicts, percentile=percentile)
        make_directory("results/motifs/filtered")
        raw_files = os.listdir("results/motifs/raw")
        for file in raw_files:
            if ".csv" in file:
                df = pd.read_csv("results/motifs/raw/%s" % file, sep=";", header=None)
                df.loc[df[4]>=self.hit_thresh].to_csv("results/motifs/filtered/%s"%file,sep=";",header=None,index=None)

    def motif_to_bed(self):
        """

        :return:
        """
        lookup_names = dict(zip(self.motif_table["Motif"], self.motif_table["Cluster_ID"]))
        motifcsvs = os.listdir("results/motifs/filtered")
        make_directory("results/motifs/bed")
        for motifcsv in motifcsvs:
            if ".csv" in motifcsv:
                df = pd.read_csv("results/motifs/filtered/%s" % motifcsv, sep=";", header=None)
                seq_len = [len(string) for string in df[5]]
                locs = df[0]
                names = df[1].values
                names = [name.split(".pmf")[0].split(" ")[0] for name in names]
                motif_ids = []
                for name in names:
                    try:
                        motif_ids.append(lookup_names[name])
                    except KeyError:
                        motif_ids.append(-1)
                chrom, start, end = [], [], []
                for loc in locs:
                    chrm, pos = loc.split(":")
                    st, en = pos.split("-")
                    st, en = int(st), int(en)
                    chrom.append(chrm)
                    start.append(st)
                    end.append(en)
                start, end = np.array(start), np.array(end)
                start += df[2].values
                end = start + seq_len
                bed_df = pd.DataFrame([chrom, start, end, motif_ids]).transpose()
                bed_df = bed_df.loc[bed_df[3] != -1]
                bed_df.to_csv("results/motifs/bed/%s.bed" % (motifcsv.split(".csv")[0]), sep="\t", header=None,
                              index=None)

                ##sort bedfile
                os.system("""
sort -k2,2n -k3,3n results/motifs/bed/%s.bed -o results/motifs/bed/%s.bed 
                """% (motifcsv.split(".csv")[0],motifcsv.split(".csv")[0]))

    def motifs_to_bedgraph(self):
        bedfiles = os.listdir("results/motifs/bed")
        make_directory("results/motifs/bedgraph")
        for bedfile in bedfiles:
            eCRE_name = bedfile.split(".bed")[0]
            bedgraph_name = eCRE_name + ".bedgraph"
            bed_to_bedgraph("results/motifs/bed/%s"%bedfile,"results/motifs/bedgraph/%s"%bedgraph_name)

    def motifs_by_archetype(self,collapse=True):
        """

        :return:
        """
        motif_beds = os.listdir("results/motifs/bed")
        make_directory("results/motifs/by_archetype")
        for bed in motif_beds:
            if ".bed" in bed:
                bedname = bed.split(".bed")[0]
                make_directory("results/motifs/by_archetype/%s"%bedname)
                beddf = pd.read_csv("results/motifs/bed/%s"%bed,sep="\t",header=None)
                for archetype in range(1,287):
                    beddf.loc[beddf[3] == archetype][beddf.columns[:3]].to_csv("results/motifs/by_archetype/%s/archetype_%d.bed"%(bedname,archetype),sep="\t",header=None,index=None)
                    if collapse is True:
                        bed_file_name = "results/motifs/by_archetype/%s/archetype_%d.bed"%(bedname,archetype)
                        BedTool(bed_file_name).merge().saveas(bed_file_name)

    def collapse_all_bed(self):
        motif_beds = os.listdir("results/motifs/bed")
        for bed in motif_beds:
            if ".bed" in bed:
                bedname = bed.split(".bed")[0]
                first = True
                for archetype in range(1,287):
                    try:
                        if first is True:
                            adf = pd.read_csv("results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype),
                                              sep="\t", header=None)
                            first = False
                        else:
                            adf = pd.concat([adf, pd.read_csv(
                                "results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype), sep="\t",
                                header=None)])
                    except:
                        a = 1
                adf.to_csv("results/motifs/bed/%s" % (bed), sep="\t", header=None,
                           index=None)


    def collapse_bed_relevant_clusters(self,clusters=[1,3,5,6]):
        """
        NB: this will duplicate archetypes found in multiple clusters

        :param clusters:
        :return:
        """
        motif_beds = os.listdir("results/motifs/bed")
        make_directory("results/motifs/relevant_clusters")
        cols = (plt.cm.Set1(np.arange(len(clusters))/len(clusters))[:,:3]*255).astype(np.int64)
        cols = [str(tuple(col)).split("(")[1].split(")")[0] for col in cols]
        for bed in motif_beds:
            if ".bed" in bed:
                bedname = bed.split(".bed")[0]
                first = True
                for i, cluster in enumerate(clusters):
                    archetypes = np.loadtxt("results/expression/archetypes/archetypes_for_cluster_%d.txt" % cluster, dtype=np.int64)
                    for archetype in archetypes:
                        try:
                            if first is True:
                                adf = pd.read_csv("results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype),
                                                  sep="\t", header=None)
                                adf[3] = "A%d"%archetype
                                score = 1000
                                adf[4] = score
                                adf[5] = "."
                                adf[6] = adf[1]
                                adf[7] = adf[2]
                                adf[8] = cols[i]
                                first = False
                            else:
                                new_df = pd.read_csv(
                                    "results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype), sep="\t",
                                    header=None)
                                new_df[3] = "A%d"%archetype
                                score = 1000
                                new_df[4] = score
                                new_df[5] = "."
                                new_df[6] = new_df[1]
                                new_df[7] = new_df[2]
                                new_df[8] = cols[i]
                                adf = pd.concat([adf, new_df])

                        except:
                            a = 1
                adf.to_csv("results/motifs/relevant_clusters/%s" % (bed), sep="\t", header=None,
                           index=None)


    def motifs_by_cluster(self,make_bedgraph=True):
        """

        :return:
        """
        archetypes_by_cluster_files = os.listdir("results/expression/archetypes")
        motif_beds = os.listdir("results/motifs/bed")
        make_directory("results/motifs/by_cluster")
        for bed in motif_beds:
            if ".bed" in bed:
                bedname = bed.split(".bed")[0]
                make_directory("results/motifs/by_cluster/%s"%bedname)
                for file in archetypes_by_cluster_files:
                    cluster_id = int(file.split(".txt")[0].split("_")[-1])
                    archetypes = np.loadtxt("results/expression/archetypes/%s"%file,dtype=np.int64)
                    first = True
                    for archetype in archetypes:
                        try:
                            if first is True:
                                adf = pd.read_csv("results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype),sep="\t",header=None)
                                first = False
                            else:
                                adf = pd.concat([adf,pd.read_csv("results/motifs/by_archetype/%s/archetype_%d.bed" % (bedname, archetype),sep="\t",header=None)])
                        except:
                            a = 1
                    adf.to_csv("results/motifs/by_cluster/%s/cluster_%d.bed"%(bedname,cluster_id),sep="\t",header=None,index=None)
                    os.system("""
    sort -k2,2n -k3,3n results/motifs/by_cluster/%s/cluster_%d.bed -o results/motifs/by_cluster/%s/cluster_%d.bed
                    """%(bedname,cluster_id,bedname,cluster_id))
                    if make_bedgraph:
                        bed_to_bedgraph("results/motifs/by_cluster/%s/cluster_%d.bed"%(bedname,cluster_id),"results/motifs/by_cluster/%s/cluster_%d.bedgraph"%(bedname,cluster_id))


class GenomePlot:
    """

    """
    def __init__(self,eCRE,plot_constructs=True,clean_configs=True,plot_bw=True,plot_genes=True,plot_phylo=True):
        self.eCRE = eCRE
        make_directory("results/genome_plots")
        make_directory("results/genome_plots/%s" % eCRE)
        make_directory("results/genome_plots/%s/config_files" % eCRE)
        make_directory("results/genome_plots/%s/config_files/all" % eCRE)
        make_directory("results/genome_plots/%s/config_files/by_cluster" % eCRE)
        make_directory("results/genome_plots/%s/config_files/by_cluster_merge" % eCRE)
        make_directory("results/genome_plots/%s/config_files/by_candidate" % eCRE)
        make_directory("results/genome_plots/%s/plots" % eCRE)
        make_directory("results/genome_plots/%s/plots/all" % eCRE)
        make_directory("results/genome_plots/%s/plots/by_cluster" % eCRE)
        make_directory("results/genome_plots/%s/plots/by_cluster_merge" % eCRE)
        make_directory("results/genome_plots/%s/plots/by_candidate" % eCRE)

        if clean_configs is True:
            self.clean_configs()

        if "_p" in self.eCRE:
            ##account for two potential eCRE inputs -- Gene_p=0.0xx and Gene
            self.eCRE_name = self.eCRE.split("_p")[0]
        else:
            self.eCRE_name = self.eCRE

        self.e_chrom, self.e_start, self.e_end = BedTool("reference/eCRE_locs/%s.bed" % self.eCRE_name).to_dataframe().values.ravel()

        self.bigwigs = pd.read_csv("reference/bigwig_files.txt", sep="\t", header=None)
        self.bigwigs.columns = ["name", "dir"]
        self.atac = pd.read_csv("reference/atac_files.txt", sep="\t", header=None)
        self.atac.columns = ["name", "dir"]
        self.phylo_files = pd.read_csv("reference/phylo_files.txt",sep="\t",header=None)
        self.phylo_files.columns = ["name", "dir"]
        self.plot_constructs = plot_constructs
        self.plot_genes = plot_genes
        self.plot_bw = plot_bw
        self.plot_phylo = plot_phylo
        if plot_constructs is True:
            self.bed_files = pd.read_csv("reference/bed_files.txt",sep="\t",header=None)
            self.bed_files.columns = ["name", "dir"]

        self.lookup = open_dict("reference/lookup_table/lookup_table")

        self.bigwig_template = """
[%s]
file=%s
title=%s
color = %s
negative_color=%s
min_value = %s
#max_value = auto
height = %.3f
number of bins = 200
nans to zeros = False
show data range = yes
file_type = bigwig
        """

        self.bigwig_template_sharey = """
        [%s]
        file=%s
        title=%s
        color = %s
        negative_color=%s
        min_value = %s
        #max_value = auto
        height = %.3f
        number of bins = 200
        nans to zeros = False
        show data range = yes
        file_type = bigwig
        overlay_previous = share-y
                """

        self.archetype_template_rows = """
[%s]
file=%s
title=%s
color = %s
height = %.3f
# line_width = 0.5
gene_rows = %d
labels = %s
file_type = bed
fontsize = 10
style = UCSC
        """

        self.archetype_template_height = """
[%s]
file=%s
title=%s
color = %s
height = %.3f
# line_width = 0.5
# gene_rows = 2
labels = %s
file_type = bed
fontsize = 10
style = UCSC
        """



        self.genes = """
[spacer]
[genes]
file = /camp/lab/luscomben/reference/Genomics/iGenomes/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf
# title of track (plotted on the right side)
title = genes
# height of track in cm (ignored if the track is overlay on top the previous track)
height = 1
# if you want to plot the track upside-down:
# orientation = inverted
# if you want to plot the track on top of the previous track. Options are 'yes' or 'share-y'.
# For the 'share-y' option the y axis values is shared between this plot and the overlay plot.
# Otherwise, each plot use its own scale
#overlay_previous = yes

# By default the transcript_name is used.
# If you want to use the gene_name:
# prefered_name = gene_name
# By default, the gtf is transformed to transcripts
# If you want to use see only one structure per gene
# merge_transcripts = true
# You can change the color of coding sequences by:
color = darkblue
# height of track in cm
height = 1
# whether printing the labels
labels = true
# optional:
# by default the labels are not printed if you have more than 60 features.
# to change it, just increase the value:
#max_labels = 60
# optional: font size can be given to override the default size
fontsize = 10
# optional: line_width
#line_width = 0.5
# the display parameter defines how the gtf file is plotted.
# Default is 'stacked' where regions are plotted on different lines so
# we can see all regions and all labels.
# The other options are ['collapsed', 'interleaved', 'triangles']
# These options assume that the regions do not overlap.
# `collapsed`: The gtf regions are plotted one after the other in one line.
# `interleaved`: The gtf regions are plotted in two lines, first up, then down, then up etc.
# optional, default is black. To remove the border, simply set 'border_color' to none
# Not used in tssarrow style
#border_color = black
# style to plot the genes when the display is not triangles
style = UCSC
#style = flybase
#style = tssarrow
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size
# otherwise, on each image the height of each gene changes
#gene_rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global_max_row = true
# If you want to plot all labels inside the plotting region:
#all_labels_inside = true
# If you want to display the name of the gene which goes over the plotted
# region in the right margin put:
#labels_in_margin = true
# if you use UCSC style, you can set the relative distance between 2 arrows on introns
# default is 2
#arrow_interval = 2
# if you use tssarrow style, you can choose the length of the arrow in bp
# (default is 4% of the plotted region)
#arrow_length = 5000
# if you use flybase or tssarrow style, you can choose the color of non-coding intervals:
#color_utr = grey
# as well as the proportion between their height and the one of coding
# (by default they are the same height):
#height_utr = 1
# By default, for oriented intervals in flybase style,
# or bed files with less than 12 columns, the arrowhead is added
# outside of the interval.
# If you want that the tip of the arrow correspond to
# the extremity of the interval use:
# arrowhead_included = true
# optional. If not given is guessed from the file ending.
file_type = gtf
        """

        self.bedgraph_template = """
[%s]
file = %s
# title of track (plotted on the right side)
title = %s
# height of track in cm (ignored if the track is overlay on top the previous track)
height = 1.5
# if you want to plot the track upside-down:
# orientation = inverted
# if you want to plot the track on top of the previous track. Options are 'yes' or 'share-y'.
# For the 'share-y' option the y axis values is shared between this plot and the overlay plot.
# Otherwise, each plot use its own scale
#overlay_previous = yes

color = green
# To use a different color for negative values
negative_color = red
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
# min_value = 0
#max_value = auto
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans_to_zeros = true
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# If you want to plot a 4C track where you want to link
# the non-missing data (NaNs) together and only use the
# middle of the region instead of the region itself:
# Default is false.
# use_middle = true
# By default the bedgraph is plotted at the base pair
# resolution. This can lead to very large pdf/svg files
# If plotting large regions.
# If you want to decrase the size of your file.
# You can either rasterize the bedgraph profile by using:
# rasterize = true
# Or use a summary method on a given number of bin:
# The possible summary methods are given by pyBigWig:
# mean/average/stdev/dev/max/min/cov/coverage/sum
# summary_method = mean
# number_of_bins = 700
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# to compute operations on the fly on the file
# or between 2 bedgraph files
# operation will be evaluated, it should contains file or
# file and second_file,
# we advice to use nans_to_zeros = true to avoid unexpected nan values
#operation = 0.89 * file
#operation = - file
#operation = file - second_file
#operation = log2((1 + file) / (1 + second_file))
#operation = max(file, second_file)
#second_file = path for the second file
# To log transform your data you can also use transform and log_pseudocount:
# For the transform values:
# 'log1p': transformed_values = log(1 + initial_values)
# 'log': transformed_values = log(log_pseudocount + initial_values)
# 'log2': transformed_values = log2(log_pseudocount + initial_values)
# 'log10': transformed_values = log10(log_pseudocount + initial_values)
# '-log': transformed_values = log(log_pseudocount + initial_values)
# For example:
#tranform = log
#log_pseudocount = 2
# When a transformation is applied, by default the y axis
# gives the transformed values, if you prefer to see
# the original values:
#y_axis_values = original
# If you want to have a grid on the y-axis
#grid = true
file_type = bedgraph
        """

        self.foot = """
[x-axis]
[spacer]
        """

        self.runline_template = """
pyGenomeTracks --tracks %s --region %s:%d-%d -o %s --width %.2f
        """

        self.runline_template_suppress = """
pyGenomeTracks --tracks %s --region %s:%d-%d -o %s --width %.2f >/dev/null 2>&1
        """

    def make_bigwig(self,name, dir, color="#666", negative_color="red",height=1.5,min_value=0,share_y=False):
        """

        :param name:
        :param dir:
        :param color:
        :param height:
        :return:
        """
        if share_y is False:
            return self.bigwig_template % (name, dir, name, color, negative_color,min_value,height)
        else:
            return self.bigwig_template_sharey % (name, dir, name, color, negative_color,min_value,height)



    def make_bed(self,name, dir, color="darkblue", height=0.75, gene_rows=None, labels="off"):
        """

        :param name:
        :param dir:
        :param color:
        :param height:
        :param gene_rows:
        :param labels:
        :return:
        """
        if gene_rows is None:
            out = self.archetype_template_height % (name, dir, name, color, height, labels)
        else:
            out = self.archetype_template_rows % (name, dir, name, color, height, gene_rows, labels)
        return out

    def write_bw(self,f,source_file=None,color="#666",min_value = 0,share_y=False):
        """

        :param f:
        :return:
        """
        if source_file is None:
            source_file = self.bigwigs
        for bwname, bwdir in source_file.values:
            if "#" not in bwname:
                f.write(self.make_bigwig(bwname, bwdir,color=color,min_value=min_value,share_y=share_y))

    def write_atac(self,f,source_file=None,colors=None,min_value = 0,share_y=False):
        """

        :param f:
        :return:
        """
        print("write atac")
        if source_file is None:
            source_file = self.atac
        if colors is None:
            colors = plt.cm.plasma(np.linspace(0,1,source_file.shape[0]))
        k = 0
        for bwname, bwdir in source_file.values:
            if "#" not in bwname:
                f.write(self.make_bigwig(bwname, bwdir,color=colors[k],min_value=min_value,share_y=share_y))
                k+=1

    def write_bedgraph(self,f):
        """

        :param f:
        :return:
        """
        for bgname,bgdir in self.bedgraph_files.values:
            if "#" not in bgname:
                f.write(self.bedgraph_template%(bgname,bgdir,bgname))

    def write_bd(self,f):
        """

        :param f:
        :return:
        """
        for bdname,bddir in self.bed_files.values:
            if "#" not in bdname:
                f.write("""
[spacer]
                """)
                f.write(self.make_bed(bdname, bddir,color="black",labels="on",height=0.3))
                f.write("""
[spacer]
                """)

    def ini_atac(self):
        """

        :return:
        """
        print("ini_atac")
        f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, "all", "atac"), 'w')
        if self.plot_bw:
            self.write_atac(f,share_y=True)
        if self.plot_genes:
            f.write(self.genes)
        if self.plot_constructs:
            self.write_bd(f)
        if self.plot_phylo:
            self.write_bw(f,self.phylo_files,color="green",min_value="auto",share_y=True)
            # self.write_bedgraph(f)
        f.write(self.foot)
        f.close()  # you can omit in most cases as the destructor will call it


    def ini_all_motifs(self):
        """

        :return:
        """
        f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, "all", "all"), 'w')
        if self.plot_bw:
            self.write_bw(f)
        if self.plot_genes:
            f.write(self.genes)
        if self.plot_constructs:
            self.write_bd(f)
        if self.plot_phylo:
            self.write_bw(f,self.phylo_files,color="green",min_value="auto")
            # self.write_bedgraph(f)
        if os.path.exists("results/motifs/bedgraph/%s.bedgraph"%(self.eCRE)):
            f.write(self.bedgraph_template%("All archetypes","results/motifs/bedgraph/%s.bedgraph"%(self.eCRE),"All archetypes"))
        else:
            f.write(self.make_bed("All archetypes", "results/motifs/bed/%s.bed" % (self.eCRE),
                             height=3))
        f.write(self.foot)
        f.close()  # you can omit in most cases as the destructor will call it

    def ini_relevant_clusters(self):
        """

        :return:
        """
        f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, "all", "relevant_clusters"), 'w')
        if self.plot_bw:
            self.write_bw(f)
        if self.plot_genes:
            f.write(self.genes)
        if self.plot_constructs:
            self.write_bd(f)
        if self.plot_phylo:
            self.write_bw(f,self.phylo_files,color="green",min_value="auto")
            # self.write_bedgraph(f)
        f.write(self.make_bed("Archetypes for relevant clusters", "results/motifs/relevant_clusters/%s.bed" % (self.eCRE),
                         height=3,color="bed_rgb"))
        f.write(self.foot)
        f.close()  # you can omit in most cases as the destructor will call it


    def ini_by_cluster_merge(self):
        """

        :return:
        """
        cat = "by_cluster_merge"
        archetype_files = os.listdir("results/motifs/by_cluster/%s" % self.eCRE)
        for archetype_file in archetype_files:
            cluster_no = int((archetype_file.split(".bed")[0]).split("cluster_")[1])
            cluster_name = "cluster_%d" % cluster_no
            f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, cat, cluster_name), 'w')
            if self.plot_bw:
                self.write_bw(f)
            if self.plot_genes:
                f.write(self.genes)
            if self.plot_constructs:
                self.write_bd(f)
            if self.plot_phylo:
                self.write_bw(f, self.phylo_files, color="green", min_value="auto")
            bedgraph_name = ("results/motifs/by_cluster/%s/%s"%(self.eCRE,archetype_file)).split(".bed")[0] + ".bedfile"
            if os.path.exists("results/motifs/by_cluster/%s/%s" % (self.eCRE,bedgraph_name)):
                f.write(self.bedgraph_template % (
                cluster_name, "results/motifs/by_cluster/%s/%s" % (self.eCRE,bedgraph_name), cluster_name))
            else:
                f.write(self.make_bed(name=cluster_name, dir="results/motifs/by_cluster/%s/%s" % (self.eCRE, archetype_file),
                                 height=3))
            f.write(self.foot)
            f.close()

    def ini_by_cluster(self):
        """

        :return:
        """
        cat = "by_cluster"
        archetype_files = os.listdir("results/motifs/by_cluster/%s" % self.eCRE)
        for archetype_file in archetype_files:
            cluster_no = int((archetype_file.split(".bed")[0]).split("cluster_")[1])
            cluster_name = "cluster_%d" % cluster_no
            f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, cat, cluster_name), 'w')
            if self.plot_bw:
                self.write_bw(f)
            if self.plot_genes:
                f.write(self.genes)
            if self.plot_constructs:
                self.write_bd(f)
            if self.plot_phylo:
                self.write_bw(f, self.phylo_files, color="green", min_value="auto")
            archetype_ids = np.loadtxt("results/expression/archetypes/archetypes_for_cluster_%d.txt" % cluster_no,
                                       dtype=np.int64)
            for aid in archetype_ids:
                f.write(self.make_bed(name="A%d" % aid,
                                 dir="results/motifs/by_archetype/%s/archetype_%d.bed" % (self.eCRE, aid),
                                 gene_rows=2))
            f.write(self.foot)
            f.close()



    def ini_by_candidate(self,
                         candidate_genes=["Nkx2-2","Nkx6-1","Irx3","Pax6","Olig2","Sox2","Gli3"]):
        """

        :param candidate_genes:
        :return:
        """
        cat = "by_candidate"
        f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (self.eCRE, cat, cat), 'w')
        if self.plot_bw:
            self.write_bw(f)
        if self.plot_genes:
            f.write(self.genes)
        if self.plot_constructs:
            self.write_bd(f)
        if self.plot_phylo:
            self.write_bw(f,self.phylo_files,color="green",min_value="auto")
        archetype_ids = [self.lookup[gene] for gene in candidate_genes]
        for i, aid in enumerate(archetype_ids):
            f.write(self.make_bed(name="%s (A%d)" % (candidate_genes[i], aid),
                             dir="results/motifs/by_archetype/%s/archetype_%d.bed" % (self.eCRE, aid),
                             gene_rows=2))
        f.write(self.foot)
        f.close()

    def make_runline(self,config_path,plot_path,suppress=True,width=20):
        """

        :param config_path:
        :param plot_path:
        :param suppress:
        :return:
        """
        if suppress is True:
            return self.runline_template_suppress%(config_path,self.e_chrom,self.e_start,self.e_end,plot_path,width)
        else:
            return self.runline_template%(config_path,self.e_chrom,self.e_start,self.e_end,plot_path,width)


    def clean_configs(self):
        for path, subdirs, files in os.walk("results/genome_plots"):
            for name in files:
                config_path = os.path.join(path, name)
                if ".ini" in config_path:
                    os.remove(config_path)

    def make_plots(self,parallel = False,suppress = True,width=20):
        """
        Makes plots for all avaliable config files
        :return:
        """

        script = """
eval "$(conda shell.bash hook)" \n
source activate pygenometracks \n
    
        """
        if parallel is True:
            scripts = []
        for path, subdirs, files in os.walk("results/genome_plots/%s/config_files"%self.eCRE):
            for name in files:
                config_path = os.path.join(path, name)
                if ".ini" in config_path:
                    plot_path = config_path.replace("config_files","plots").replace(".ini",".pdf")
                    runline = self.make_runline(config_path,plot_path,suppress=suppress,width=width)
                    if parallel is False:
                        script += runline
                    else:
                        scripts.append(script+runline)

        if parallel is False:
            os.system(script)
        else:
            num_cores = multiprocessing.cpu_count()
            Parallel(n_jobs=num_cores)(delayed(os.system)(scriptt) for scriptt in scripts)

def make_directory(dir):
    """

    :param dir:
    :return:
    """
    if not os.path.exists(dir):
        os.mkdir(dir)

def delete_directory(dir):
    """

    :param dir:
    :return:
    """
    if os.path.exists(dir):
        shutil.rmtree(dir)

def save_csv(data,filename,**kwargs):
    """

    :param data:
    :param filename:
    :param kwargs:
    :return:
    """
    gene_df = pd.DataFrame(data)
    gene_df.to_csv(filename,**kwargs)

def capitalize_list(list):
    """

    :param list:
    :return:
    """
    return [val.capitalize() for val in list]

def save_dict(dictionary,filename):
    """

    :param dictionary:
    :param filename:
    :return:
    """
    with open('%s.json'%filename, 'w') as fp:
        json.dump(dictionary, fp)

def open_dict(filename):
    """

    :param filename:
    :return:
    """
    with open('%s.json'%filename, 'r') as fp:
        dictionary = json.load(fp)
    return dictionary

def standardize_name(name,aliases,true_names=None):
    """

    :param name:
    :param aliases:
    :param true_names:
    :return:
    """
    if true_names is None:
        true_names = aliases.values()
    if name not in true_names:
        try:
            return aliases[name]
        except KeyError:
            return "*" + name
    else:
        return name

def standardize_names(names,aliases,true_names=None):
    """

    :param names:
    :param aliases:
    :param true_names:
    :return:
    """
    return [standardize_name(name, aliases,true_names) for name in names]


def unzip(file,outfile):
    """

    :param file:
    :param outfile:
    :return:
    """
    with gzip.open(file, 'r') as f_in, open(outfile, 'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)

def bed_to_bedgraph(input,output,genome="mm10"):
    """

    :return:
    """
    BedTool(input).genomecov(bg=True,genome=genome).saveas(output)
    print("Converted",input,"to",output)