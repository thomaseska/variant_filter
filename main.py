import functions.parse_vcf as parse_vcf
import functions.call_GN as call_GN
import functions.parse_gn_result as parse_gn_result
import functions.helpers as helpers
import argparse
import json
import os

def run(proxies, verify,
        vcf_path, excel_file, json_out,
        min_qual, min_dp, min_VF, max_gnomad,
        fields, token):


    req_ls, var_freq = parse_vcf.parse_vcf(vcf_path=vcf_path,
                                            min_qual=min_qual,
                                            min_dp=min_dp,
                                            min_VF=min_VF)

    gnomad = call_GN.request_variant_info(req_ls=req_ls,
                                        proxies=proxies,
                                        verify=verify,
                                        max_gnomad=max_gnomad)

    req_ls = [variant for variant in req_ls if variant in list(gnomad.keys())]

    print(f"{helpers.nice_time()} : {len(req_ls)} variants left after gnomad filtering")

    resp_ls_dics = call_GN.request_annotation(req_ls=req_ls,
                                            proxies=proxies,
                                            verify=verify,
                                            fields= fields,
                                            token=token,
                                            json_out=json_out)



    out_df = parse_gn_result.parse_result(resp_ls_dics, var_freq, gnomad)

    out_df.to_excel(excel_file, index=False)

    print(f"{helpers.nice_time()} : Run complete")



parser = argparse.ArgumentParser(description='Process VCF files and retrieve variant annotation')

parser.add_argument('VCF_file_or_folder',
                     help="VCF file to process. If running with -b, folder contatining VCF files")
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


min_qual = args.quality
min_dp = args.depth
min_VF = args.frequency
max_gnomad = args.gnomad

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
    
    run(vcf_path=vcf_path, excel_file=excel_file, json_out=json_out,
        proxies=proxies, verify=verify, token=token, fields=fields,
        min_qual=min_qual, min_dp=min_dp, min_VF=min_VF, max_gnomad=max_gnomad)

else:
    vcf_folder = args.VCF_file_or_folder
    out_folder = args.output
    vcf_files = os.listdir(vcf_folder)

 
    for file in vcf_files:
        vcf_path = f"{vcf_folder}/{file}"
        excel_file = f"{out_folder}/{file[:-4]}.xlsx"

        if args.json:
            json_out = f"{args.json}/{file[:-4]}.json"
        else:
            json_out = None

        run(vcf_path=vcf_path, excel_file=excel_file, json_out=json_out,
            proxies=proxies, verify=verify, token=token, fields=fields,
            min_qual=min_qual, min_dp=min_dp, min_VF=min_VF, max_gnomad=max_gnomad)
        continue
        # start runs

# vcf_path = "files/E1130_var_only.vcf"
# vcf_path = "files/E1130_23Pabcd_MergedSmallVariants.genome.vcf"
# vcf_path = "files/test.vcf"

