
configfile: "config.yaml"

INPUT=glob_wildcards(f"{config['base_dir']}/{{file}}.vcf").file

print(INPUT)

rule all:
    input:
        expand(f"{config['result_dir']}/{{name}}_q{config['quality']}.xlsx", name=INPUT)


rule run_pipeline:
    input:
        f"{config['base_dir']}/{{sample}}.vcf"
    output:
        xlsx=f"{config['result_dir']}/{{sample}}_q{config['quality']}.xlsx",
        csv=f"{config['result_dir']}/{{sample}}_q{config['quality']}.csv"
    log:
        f"{config['result_dir']}/{{sample}}_q{config['quality']}.log"

    run:
        if config["rmflag"]:
            if config["mutect2"]:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    --mutect2 \
                    --rmflag {config[rmflag]} \
                    --rmflagfile /sound/data/genomes/hsapiens/flag_genes_list.txt \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ")
            elif config["ont"]:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    --ont \
                    --rmflag {config[rmflag]} \
                    --rmflagfile /sound/data/genomes/hsapiens/flag_genes_list.txt \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ")
            else:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    --rmflag {config[rmflag]} \
                    --rmflagfile /sound/data/genomes/hsapiens/flag_genes_list.txt \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ") 
        else:
            if config["mutect2"]:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    --mutect2 \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ")
            elif config["ont"]:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    --ont \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ")
            else:
                shell("python /sound/home/thomas/variant_filter/main.py \
                    -o {output.xlsx} \
                    --csv {output.csv} \
                    --proxies /sound/home/thomas/variant_filter/proxies.json \
                    --ssl /etc/ssl/certs \
                    --token /sound/home/thomas/variant_filter/oncoKB_token.json \
                    -d {config[depth]} \
                    -q {config[quality]} \
                    -f {config[frequency]} \
                    -g {config[gnomad]} \
                    {input} &> {log} \
                ") 

