INPUT=glob_wildcards("vcf_only/vcfs_hg19/{file}.vcf").file

rule all:
    input:
        expand("filtered_variant_files/{name}.xlsx", name=INPUT)


rule run_pipeline:
    input:
        "vcf_only/vcfs_hg19/{sample}.vcf"
    output:
        xlsx="filtered_variant_files/{sample}.xlsx",
        csv="filtered_variant_files/{sample}.csv"
    log:
        "filtered_variant_files/{sample}.log"
    shell:
        "python /sound/home/thomas/variant_filter/main.py -o {output.xlsx} --csv {output.csv} --proxies /sound/home/thomas/variant_filter/proxies.json --ssl /etc/ssl/certs"
        " --token /sound/home/thomas/variant_filter/oncoKB_token.json --mutect2 --rmflag 100 --rmflagfile /sound/data/genomes/hsapiens/flag_genes_list.txt -d 2"
        " {input} &> {log}"
