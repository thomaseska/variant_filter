import pandas as pd
from . import helpers

def parse_patho(csv_path):
    if csv_path[-4:] == ".csv":
        df = pd.read_csv(csv_path, sep=";")
        if len(df.columns) == 1:
            df = pd.read_csv(csv_path, sep=",")
            if len(df.columns) == 1:
                print("Input is not separated by ',' or ';'. Please check your input file!")
    elif csv_path[-5:] == ".xlsx":
        df = pd.read_excel(csv_path)
    else:
        print("Input is not .csv or .xlsx. Please check your input file!")
        exit(1)
    
    # select case:
    if "pathoId" in df.columns.values.tolist():
        unique_cases = df["pathoId"].unique()
        print("Available IDs:")
        for i, id in enumerate(unique_cases, start=1):
            print(f"[{i}] {id}")
        choice = int(input("Enter your Selection: "))

        selected_value = unique_cases[choice -1]
        df = df[df["pathoId"] == selected_value]
    elif "ng_id" in df.columns.values.tolist():
        unique_cases = df["ngs_id"].unique()
        print("Available IDs:")
        for i, id in enumerate(unique_cases, start=1):
            print(f"[{i}] {id}")
        choice = int(input("Enter your Selection: "))

        selected_value = unique_cases[choice -1]
        df = df[df["ngs_id"] == selected_value]
    else:
        print("Identifier Column not found. Select Identifier column from options below:")
        print("[0] No Identifier (Use all rows)")
        for i, id in enumerate(df.columns.values.tolist(), start=1):
            print(f"[{i}] {id}")
        col_choice = int(input("Enter your Selection: "))
        if col_choice != 0:
            col_selected_value = df.columns.values.tolist()[col_choice -1]
            
            unique_cases = df[col_selected_value].unique()
            print("Available IDs:")
            for i, id in enumerate(unique_cases, start=1):
                print(f"[{i}] {id}")
            choice = int(input("Enter your Selection: "))

            selected_value = unique_cases[choice -1]
            df = df[df[col_selected_value] == selected_value]

    print(df)

    df = df[["chr", "start", "ref", "alt"]]
    df = df.dropna()
    print(df)
    def build_hgvsg(row):
        if len(row["ref"]) > 1:
            if len(row["alt"]) > 1:
                # deletion-insertion
                hgvsg = f'{row["chr"][3:]}:g.{int(row["start"])}_{int(row["start"]) + len(row["alt"]) - 1}delins{row["alt"]}'
            # deletion
            else:
                if len(row["ref"]) > 2:
                    hgvsg = f'{row["chr"][3:]}:g.{int(row["start"]) + 1}_{int(row["start"]) + len(row["ref"]) - 1}del'
                else:
                    hgvsg = f'{row["chr"][3:]}:g.{int(row["start"]) + 1}del'
        # insertion
        elif len(row["alt"]) > 1:
                hgvsg = f'{row["chr"][3:]}:g.{int(row["start"])}_{int(row["start"]) + 1}ins{row["alt"][1:]}'
        # substitution
        else:
            hgvsg = f'{row["chr"][3:]}:g.{row["start"]}{row["ref"]}>{row["alt"]}'
        return hgvsg
    
    hgvsg_ls = list(df.apply(build_hgvsg, axis=1))
    
    return hgvsg_ls



def parse_patho_result(resp_ls_dics):
    """
    Filter genome nexus results and retain (likely) clinically relevant variant transcripts.
    Filters are ClinVar, oncoKB, PolyPhen and SIFT terms.
    """
    print(f"{helpers.nice_time()} : Creating excel file")

    out_dict = []
    for resp_dic in resp_ls_dics:
        try:
            hgvsg = resp_dic["hgvsg"]
        except KeyError:
            continue
        # clinvar terms
        try:
            clinvar = resp_dic["clinvar"]["annotation"]["clinicalSignificance"]
        except KeyError:
            clinvar = ""
        except TypeError:
            clinvar = ""

        # oncokb
        try:
            oncokb1 = resp_dic["oncokb"]["annotation"]["mutationEffect"]["knownEffect"]
        except KeyError:
            oncokb1 = ""
        try:
            oncokb2 = resp_dic["oncokb"]["annotation"]["oncogenic"]
        except KeyError:
            oncokb2 = ""
        
        if oncokb1 == "":
            oncokb = oncokb2
        elif oncokb2 == "":
            oncokb == oncokb1
        else:
            oncokb = f"{oncokb1},{oncokb2}"

        # if there are no transcripts and terms are tolerated return the variant
        if len(resp_dic["annotation_summary"]["transcriptConsequenceSummaries"]) == 0:
 
            dic_to_add = {
                "hgvsg": hgvsg,
                "clinvar": clinvar,
                "oncokb": oncokb
            }

            out_dict.append(dic_to_add)            
        
        # loop through transcripts
        for conseq in resp_dic["annotation_summary"]["transcriptConsequenceSummaries"]:
            
            transcript_id = conseq["transcriptId"]
            hugo = conseq["hugoGeneSymbol"]
            if "hgvsc" in conseq:
                hgvsc = conseq["hgvsc"]
            else:
                hgvsc = ""
            
            if "hgvsp" in conseq:
                hgvsp = conseq["hgvsp"]
            else:
                hgvsp = ""

            conseq_terms = conseq["consequenceTerms"]
            variant_class = conseq["variantClassification"]

            # polyphen
            if "polyphenScore" in conseq:
                poly_score = conseq["polyphenScore"]
                poly_pred = conseq["polyphenPrediction"]
            else:
                poly_score = ""
                poly_pred = ""
            
            # sift
            if "siftScore" in conseq:
                sift_score = conseq["siftScore"]
                sift_pred = conseq["siftPrediction"]
            else:
                sift_score = ""
                sift_pred = ""

            
            # add transcript to output
            dic_to_add = {
                "hugo_gene_symbol": hugo,
                "hgvsg": hgvsg,
                "hgvsc": hgvsc,
                "hgvsp": hgvsp,
                "consequence_terms": conseq_terms,
                "variant_classification": variant_class,
                "clinvar": clinvar,
                "oncokb": oncokb,
                # tmp
                "sift_score": sift_score,
                "sift_prediction": sift_pred,
                "polyphen_score": poly_score,
                "polyphen_prediction": poly_pred
            }

            out_dict.append(dic_to_add)

    df = pd.DataFrame.from_records(out_dict)
    df = df.drop_duplicates()
    return df


