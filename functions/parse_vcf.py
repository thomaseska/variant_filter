from . import helpers

def parse_vcf(vcf_path, min_qual=50, min_dp=25, min_VF=0.05):
    """
    Receive a VCF formatted file and return variant HGVSg strings and variant frequency.
    Filters:    min_qual    | minimum required quality score
                min_dp      | minimum required sequencing depth
                min_VF      | minimum required variant frequency within the sample
    Created for files generated by TSO500.

    Mutations on Blacklist will be removed. Mutations in TERT and a short region downstream will be retained, regardless of depth and quality.
    """

    print(f"Assessing file {vcf_path}")
    print(f"{helpers.nice_time()} : Begin processing...")

    linecounter = 0
    format = None

    req_ls = []
    var_freq = {}

    with open(vcf_path, "r") as file:
        for line in file:
            linecounter += 1
            if linecounter % 100000 == 0:
                print(f"{helpers.nice_time()} : Processed {linecounter} lines")

            # filter comments
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            # filter non variants
            if fields[4] == ".":
                continue

            # filter quality
            if int(fields[5]) < min_qual:
                # keep TERT
                if not (fields[0] == "chr5" and int(fields[1]) >= 1253147 and int(fields[1]) <= 1295300):
                    continue

            # filter depth
            if int(fields[7].split("=")[1]) < min_dp:
                # keep TERT
                if not (fields[0] == "chr5" and int(fields[1]) >= 1253147 and int(fields[1]) <= 1295300):
                    continue

            # filter allele frequency
            format = fields[8].split(":")
            vf_loc = format.index("VF")
            values = fields[9].split(":")
            if float(values[vf_loc]) < min_VF:
                continue

            # filter blacklist
            if "Blacklist" in fields[6].split(";"):
                continue

            # build identifier
            ref = fields[3]
            alt = fields[4]
            if len(ref) > 1:
                if len(alt) > 1:
                    # deletion-insertion
                    hgvsg = f"{fields[0][3:]}:g.{int(fields[1])}_{int(fields[1]) + len(alt) - 1}delins{alt}"
                # deletion
                else:
                    if len(ref) > 2:
                        hgvsg = f"{fields[0][3:]}:g.{int(fields[1]) + 1}_{int(fields[1]) + len(ref) - 1}del"
                    else:
                        hgvsg = f"{fields[0][3:]}:g.{int(fields[1]) + 1}del"
            # insertion
            elif len(alt) > 1:
                    hgvsg = f"{fields[0][3:]}:g.{int(fields[1])}_{int(fields[1]) + 1}ins{alt[1:]}"
            # substitution
            else:
                hgvsg = f"{fields[0][3:]}:g.{fields[1]}{ref}>{alt}"

            req_ls.append(hgvsg)
            var_freq[hgvsg] = values[vf_loc]

    print(f"{helpers.nice_time()} : Found {len(req_ls)} relevant variants according to filtering parameters")
    return(req_ls, var_freq)
