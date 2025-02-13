from gwf import Workflow, AnonymousTarget
# from gwf.workflow import collect
import pandas as pd
import os
from groups import Group

gwf = Workflow()

# Inputs for running Relate, both with autosomal data, and female_only X
path_to_relate = "/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/relate/relate_v1.1.7_x86_64_dynamic/"
path_to_vcfs = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.phased.rehead" #  Autosomal
path_to_vcfs_females = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.phased.females" #  Females only
path_to_ancestor = "/home/eriks/baboondiversity/data/ancestral_state_panu3_23_04_2021/papio_anubis_ancestor_{}.fa"
path_to_poplabels = "data/pops/all_inds_8cluster.sample"  # Autosomal
path_to_poplabels_females = "data/pops/all_females_8cluster.sample"  # Females only, for X
path_to_mask = "/home/eriks/baboondiversity/data/callability_panu3_26_04_2021/panu3.npmask.chr{}.fa"
genetic_map = "/home/eriks/baboondiversity/data/PG_panu3_recombination_map/mikumi_pyrho_genetic_map_chr{}.txt"
haps_sample_dir = "steps/haps_sample/"

pop_files = "data/pops/"
keep_pop_l_l = [["Eastern Yellow"], ["All_Samples"]]
pop_information = ["0.57e-8", "32000"]

vcfs = []
autosomes = list(range(19, 21))
chromosomes = autosomes + ["X"]

for chrom in autosomes:
    vcf_path_and_name = os.path.join(path_to_vcfs.format(chrom, chrom))
    vcfs.append({"vcf_path": vcf_path_and_name, "chrom": chrom})
x_vcf = [{"vcf_path": path_to_vcfs_females.format("X", "X"), "chrom": "X"}]

l_d = []
for chrom in chromosomes:
    d = {}
    d["number"] = chrom
    d["genetic_map"] = genetic_map.format(chrom)
    l_d.append(d)


def vcf_to_haps(vcf_path, chrom, relate_path, output_dir):
    """Converts vcf files to haps/sample using the script from Relate """
    inputs = [vcf_path+".vcf.gz"]
    haps_out = os.path.join(output_dir, "chrom{}.haps". format(chrom))
    sample_out = os.path.join(output_dir, "chrom{}.sample".format(chrom))
    RelateFileFormats = os.path.join(relate_path, "bin/RelateFileFormats")
    outputs = {"haps": haps_out, "sample": sample_out}
    options = {
        "cores": 2,
        "memory": "4g",
        "walltime": "1:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} --mode ConvertFromVcf --haps {} --sample {} -i {}
    """.format(RelateFileFormats, haps_out, sample_out, vcf_path)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prepare_input(haps, sample, mask, ancestor, pop_l, output_dir, poplabels, relate_path):
    inputs = {"haps": haps, "sample": sample}
    s = os.path.basename(haps)
    number = s[5:s.find(".")]
    PrepareInputFiles = relate_path + "scripts/PrepareInputFiles/PrepareInputFiles.sh"
    destination_name = output_dir + "chrom{}".format(number)
    n_mask = mask.format(number)
    ancestor_in = ancestor.format(number)
    outputs = {"haps": destination_name+".haps.gz", "sample": destination_name+".sample.gz"}
    remove_ids = ""
    pop_file = pd.read_csv(poplabels, sep=" ")
    if pop_l != ["All_Samples"]:
        ids_to_remove = pop_file.loc[~pop_file["POP"].isin(pop_l)]
        ids_to_remove.ID.to_csv(output_dir+"IDs_removed{}.txt".format(number),index=False, header=False)
        remove_ids = " --remove_ids "+output_dir+"IDs_removed{}.txt".format(number)
        labels_to_keep = pop_file.loc[pop_file["POP"].isin(pop_l)]
        labels_to_keep.to_csv(destination_name+".poplabels", index=False, sep=" ")
    else:
        pop_file.to_csv(destination_name+".poplabels", index=False, sep=" ")
    options = {
        "cores": 2,
        "memory": "4g",
        "walltime": "4:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {}  --haps {} --sample {} --mask {} --ancestor {}{} --poplabels {} -o {}
    """.format(PrepareInputFiles, haps, sample, n_mask, ancestor_in, remove_ids, poplabels, destination_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def full_relate(number, genetic_map, out_dir, pop_inf, prep_dir, base_path, relate_path):
    haps = prep_dir+"chrom{}.haps.gz".format(number)
    sample = prep_dir+"chrom{}.sample.gz".format(number)
    annot = prep_dir+"chrom{}.annot".format(number)
    dist = prep_dir+"chrom{}.dist.gz".format(number)
    m = pop_inf[0]
    n = pop_inf[1]
    o = "chrom{}".format(number)
    inputs = {"haps": haps, "sample": sample}
    Relate = os.path.join(relate_path, "bin/Relate")
    outputs = [out_dir+o+".anc"]
    options = {
        "cores": 4,
        "memory": "40g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    cd {}
    {} --mode All -m {} -N {} --haps {} --sample {} --map {} -o {} --annot {} --dist {} --memory 20
    """.format(out_dir,
               Relate, m, n, haps, sample, genetic_map, o, annot,  dist)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def estimate_pop_size(number, out_dir, pop_inf, prep_dir, relate_path):
    i = out_dir+"chrom{}.anc".format(number)
    o = out_dir+"chrom{}_popsize".format(number)
    m = pop_inf[0]
    inputs = [out_dir+"chrom{}.anc".format(number)]
    outputs = [o+".coal"]
    Relate = os.path.join(relate_path, "scripts/EstimatePopulationSize/EstimatePopulationSize.sh")
    poplabels = prep_dir+"chrom{}.poplabels".format(number)
    options = {
        "cores": 14,
        "memory": "60g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} -i {} -m {} --poplabels {} -o {} --threshold 0 --num_iter 5 --years_per_gen 11 --threads 14
    """.format(Relate, i[:-4], m, poplabels, o)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def detect_selection(number, in_dir, out_results, pop_inf, prep_dir, relate_path):
    i = in_dir+"chrom{}_popsize.coal".format(number)
    o = out_results+"chrom{}_selection".format(number)
    m = pop_inf[0]
    inputs = [i]
    outputs = [o+".freq"]
    Relate = os.path.join(relate_path, "scripts/DetectSelection/DetectSelection.sh")
    options = {
        "cores": 5,
        "memory": "20g",
        "walltime": "1:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} -i {} -m {} -o {}
    """.format(Relate, i[:-5], m, o)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Setting up and running the workflow - autosome and X are split due to having different
# sample sizes and metadata.

os.makedirs(haps_sample_dir, exist_ok=True)
haps_sample_aut = gwf.map(vcf_to_haps, vcfs, extra={"relate_path": path_to_relate, "output_dir": haps_sample_dir})
haps_sample_x = gwf.map(vcf_to_haps, x_vcf, extra={"relate_path": path_to_relate, "output_dir": haps_sample_dir},
                        name='x_to_haps')

for pop_l in keep_pop_l_l:
    popname = "_".join([x.replace(" ", "_") for x in pop_l])
    print(popname)
    path_to_prep = os.getcwd()+"/steps/relate_{}_prepared/".format(popname)
    os.makedirs(path_to_prep, exist_ok=True)
    os.makedirs(path_to_prep+"temp/", exist_ok=True)
    with Group(gwf, suffix=popname) as g:
        preps_aut = g.map(prepare_input, haps_sample_aut.outputs, name="prepare_input_aut", extra={
                        "mask": path_to_mask, "ancestor": path_to_ancestor, "pop_l": pop_l,
                        "output_dir": path_to_prep, "poplabels": path_to_poplabels,
                        "relate_path": path_to_relate})
        preps_x = g.map(prepare_input, haps_sample_x.outputs, name="prepare_input_x", extra={
                        "mask": path_to_mask, "ancestor": path_to_ancestor, "pop_l": pop_l,
                        "output_dir": path_to_prep, "poplabels": path_to_poplabels_females,
                        "relate_path": path_to_relate})
    out_dir = "steps/{}_relate/".format(popname)
    out_results = "results/{}_relate/".format(popname)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(out_results, exist_ok=True)
    with Group(gwf, suffix=popname) as g:
        relate = g.map(full_relate, l_d, name="relate", extra={
            "out_dir": out_dir, "pop_inf": pop_information,
            "prep_dir": path_to_prep, "base_path": os.getcwd(), "relate_path": path_to_relate
        })
        popsize = g.map(estimate_pop_size, chromosomes, name="popsize", extra={
            "out_dir": out_dir, "pop_inf": pop_information,
            "prep_dir": path_to_prep, "relate_path": path_to_relate
        })
        g.map(detect_selection, chromosomes, name="detect_selection", extra={
            "in_dir": out_dir, "out_results": out_results, "pop_inf": pop_information,
            "prep_dir": path_to_prep, "relate_path": path_to_relate
        })
