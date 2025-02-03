import functions.parse_vcf as parse_vcf
import functions.call_GN as call_GN
import functions.parse_gn_result as parse_gn_result
import functions.helpers as helpers
import argparse
import json
import os
import glob
import pandas as pd

def run(proxies, verify,
        vcf_path, excel_file, json_out, csv_out,
        min_qual, min_dp, min_VF, max_gnomad, flags,
        fields, token, ont=False, mutect=False):
    """
    Execute the workflow for one VCF file:
        - parse the VCF to obtain relevant variants
        - get population frequency of variants from myvariant.info/gnomAD
        - remove common variants
        - obtain clinical information from genome nexus
        - retain clinically relevant variant transcripts 
    """

    if ont:
        req_ls, var_freq = parse_vcf.parse_vcf_ont(vcf_path=vcf_path,
                                                    min_qual=min_qual,
                                                    min_dp=min_dp,
                                                    min_VF=min_VF,
                                                    flags=flags)
    elif mutect:
        req_ls, var_freq = parse_vcf.parse_vcf_mutect(vcf_path=vcf_path,
                                            pass_filter = ["PASS"],
                                            min_dp=min_dp,
                                            min_VF=min_VF,
                                            flags=flags)

    else:
        req_ls, var_freq = parse_vcf.parse_vcf(vcf_path=vcf_path,
                                                min_qual=min_qual,
                                                min_dp=min_dp,
                                                min_VF=min_VF,
                                                flags=flags)

    if not mutect:
        gnomad = call_GN.request_variant_info(req_ls=req_ls,
                                            proxies=proxies,
                                            verify=verify,
                                            max_gnomad=max_gnomad)

        # filter variant list by gnomAD result
        req_ls = [variant for variant in req_ls if variant in list(gnomad.keys())]
    else:
        gnomad = {}

    print(f"{helpers.nice_time()} : {len(req_ls)} variants left after gnomad filtering")

    resp_ls_dics = call_GN.request_annotation(req_ls=req_ls,
                                            proxies=proxies,
                                            verify=verify,
                                            fields= fields,
                                            token=token,
                                            json_out=json_out)



    out_df = parse_gn_result.parse_result(resp_ls_dics, var_freq, gnomad)

    # if len(flags) > 0:
    #     out_df = out_df[~out_df["hugo_gene_symbol"].isin(flags)]

    try:
        out_df.to_excel(excel_file, index=False)
    except ValueError as e:
        print(f"{helpers.nice_time()} : ERROR creating xlsx file ({e})")


    if csv_out:
        out_df.to_csv(csv_out, index=False)

    print(f"{helpers.nice_time()} : Run complete")



parser = argparse.ArgumentParser(description='Process VCF files and retrieve variant annotation')

parser.add_argument('VCF_file_or_folder',
                     help="VCF file to process. If running with -b, folder contatining only VCF files")
parser.add_argument("-o", "--output", help="Output file to write to (xlsx)", required=True)
parser.add_argument("-j", "--json", help="file to save Genome Nexus json response to")

parser.add_argument("-b", "--batch", action="store_true", help="execute in batch mode")
parser.add_argument("-p", "--proxies", help="path to JSON containing proxies")
parser.add_argument("-s", "--ssl", help="path to SSL certificate folder")

parser.add_argument("-q", "--quality",
                     help="Minimum quality score for a variant to be considered (Default 50)",
                     default=50, type=int)
parser.add_argument("-d", "--depth",
                    help="Minimum allele depth for a variant to be considered (Default 25)",
                    default=25, type=int)
parser.add_argument("-f", "--frequency",
                    help="Minimum variant frequency in the sample for the variant to be considered (Default 0.05)",
                    default=0.05, type=float)
parser.add_argument("-g", "--gnomad",
                    help="Maximum variant frequency in the population for the variant to be considered (Default 0.01)",
                    default=0.01, type=float)

parser.add_argument("-t", "--token", help="file containing token(s)")

parser.add_argument("--ont", help="VCF was generated from ONT data", action="store_true")

parser.add_argument("--mutect2", help="VCF was generated using mutect2", action="store_true")

parser.add_argument("--csv", help="optional output file in csv format")

parser.add_argument("--rmflag", help="remove the top N flag genes from the output (default 0)",
                    default=0, type=int)

# https://static-content.springer.com/esm/art%3A10.1186%2Fs12920-017-0309-7/MediaObjects/12920_2017_309_MOESM3_ESM.txt
parser.add_argument("--rmflagfile", help="path to file containing flag genes")

args = parser.parse_args()

if args.proxies:
    with open(args.proxies) as prox:
        proxies = json.load(prox)
else:
    proxies = None

if args.ssl:
    verify = args.ssl
else:
    verify = None

if args.ont:
    ont = True
else:
    ont = False

if args.mutect2:
    mutect = True
else:
    mutect = False

if mutect & ont:
    print("Options ont and mutect2 are mutually exclusive! Exiting...")
    exit(1)



min_qual = args.quality
min_dp = args.depth
min_VF = args.frequency
max_gnomad = args.gnomad
rmflag = args.rmflag


flags = None
if rmflag > 0:
    try:
        # counter = 1
        # with open(args.rmflagfile, "r")as flagfile:
        #     for line in flagfile:
        #         hugo = line.strip("\n").split("\t")[0]
        #         flags.append(hugo)
        #         if counter == rmflag:
        #             break
        flags = pd.read_csv(args.rmflagfile, names = ["hgnc", "freq", "chr", "start", "stop", "hgnc_num"], sep=",")
        flags = flags.iloc[0:rmflag,]

    except FileNotFoundError:
        print("Flagfile not found. Exiting...")
        exit(0)


fields = ["annotation_summary", "clinvar", "oncokb", "my_variant_info"]

if args.token:
    with open(args.token) as tok:
        token = json.load(tok)
else:
    token = None


print(f"Begin Variant Filtering...")
print(f"Filtering parameters: ")
print(f"    - Minimum quality score: {min_qual}")
print(f"    - Minimum sequencing depth: {min_dp}")
print(f"    - Minimum allele frequency: {min_VF}")
print(f"    - Maximum population allele frequency: {max_gnomad}")

if not args.batch:
    vcf_path = args.VCF_file_or_folder
    excel_file = args.output
    if args.json:
        json_out = args.json
    else:
        json_out = None
    
    if args.csv:
        csv_out = args.csv
    else:
        csv_out = None
    
    run(vcf_path=vcf_path, excel_file=excel_file, json_out=json_out, csv_out=csv_out,
        proxies=proxies, verify=verify, token=token, fields=fields,
        min_qual=min_qual, min_dp=min_dp, min_VF=min_VF, max_gnomad=max_gnomad, flags=flags,
        ont=ont, mutect=mutect)

else:
    vcf_folder = args.VCF_file_or_folder
    out_folder = args.output
    vcf_files = glob.glob("*.vcf", root_dir=vcf_folder)

 
    for file in vcf_files:
        vcf_path = f"{vcf_folder}/{file}"
        excel_file = f"{out_folder}/{file[:-4]}.xlsx"

        if args.json:
            json_out = f"{args.json}/{file[:-4]}.json"
        else:
            json_out = None
        
        if args.csv:
            csv_out = f"{args.csv}/{file[:-4]}.csv"
        else:
            csv_out = None

        run(vcf_path=vcf_path, excel_file=excel_file, json_out=json_out, csv_out=csv_out,
            proxies=proxies, verify=verify, token=token, fields=fields,
            min_qual=min_qual, min_dp=min_dp, min_VF=min_VF, max_gnomad=max_gnomad, flags=flags,
            ont=ont, mutect=mutect)
