from . import helpers
import pandas as pd

def parse_result(resp_ls_dics, var_freq, gnomad):
    print(f"{helpers.nice_time()} : Creating excel file")

    clinvar_tolerated = ["Benign", "Likely_benign"]
    oncokb_tolerated = [ "Neutral", "Likely Neutral", "Unknown"] # 
    polyphen_tolerated = ["unknown", "benign"]
    sift_tolerated = ["tolerated"]
    out_dict = []
    for resp_dic in resp_ls_dics:
        hgvsg = resp_dic["hgvsg"]
        try:
            clinvar = resp_dic["clinvar"]["annotation"]["clinicalSignificance"]
        except KeyError:
            clinvar = ""
        except TypeError:
            clinvar = ""

        # gnomad/allele frequency
        try:
            pop_frequency = gnomad[hgvsg]["gnomad_genome"]["af"]["af"]
        except KeyError:
            pop_frequency = ""
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

            if "polyphenScore" in conseq:
                poly_score = conseq["polyphenScore"]
                poly_pred = conseq["polyphenPrediction"]
            else:
                poly_score = ""
                poly_pred = ""

            if "siftScore" in conseq:
                sift_score = conseq["siftScore"]
                sift_pred = conseq["siftPrediction"]
            else:
                sift_score = ""
                sift_pred = ""
            
            if clinvar == "" or all(x in clinvar_tolerated for x in clinvar.split("/")):

                if oncokb == "" or all(x in oncokb_tolerated for x in oncokb.split(",")):
                    
                    if sift_score == "" or sift_pred in sift_tolerated:

                        if poly_pred == "" or poly_pred in polyphen_tolerated:

                            continue
            

            dic_to_add = {
                "hugo_gene_symbol": hugo,
                "hgvsg": hgvsg,
                "hgvsc": hgvsc,
                "hgvsp": hgvsp,
                "consequence_terms": conseq_terms,
                "variant_classification": variant_class,
                "allele_frequency": var_freq[hgvsg],
                "clinvar": clinvar,
                "oncokb": oncokb,
                # tmp
                "sift_score": sift_score,
                "sift_prediction": sift_pred,
                "polyphen_score": poly_score,
                "polyphen_prediction": poly_pred,
                "population_frequency": pop_frequency
            }

            out_dict.append(dic_to_add)

    df = pd.DataFrame.from_records(out_dict)

    return df

