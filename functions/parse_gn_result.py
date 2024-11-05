from . import helpers
import pandas as pd

def parse_result(resp_ls_dics, var_freq, gnomad):
    """
    Filter genome nexus results and retain (likely) clinically relevant variant transcripts.
    Filters are ClinVar, oncoKB, PolyPhen and SIFT terms.
    """
    print(f"{helpers.nice_time()} : Creating excel file")

    # tolerated terms for each information source
    clinvar_tolerated = ["Benign", "Likely_benign"]
    oncokb_tolerated = [ "Neutral", "Likely Neutral", "Unknown"]
    polyphen_tolerated = ["benign"]
    sift_tolerated = ["tolerated"]

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

        # if there are no transcripts and terms are tolerated return the variant
        if len(resp_dic["annotation_summary"]["transcriptConsequenceSummaries"]) == 0:
                
            if all(x in clinvar_tolerated for x in clinvar.split("/")):
                if all(x in oncokb_tolerated for x in oncokb.split(",")):
                    continue
 
            dic_to_add = {
                "hgvsg": hgvsg,
                "allele_frequency": var_freq[hgvsg],
                "clinvar": clinvar,
                "oncokb": oncokb,
                "population_frequency": pop_frequency
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
            
            # if all terms are tolerated, skip the transcript
            if clinvar == "" or all(x in clinvar_tolerated for x in clinvar.split("/")):

                if oncokb == "" or all(x in oncokb_tolerated for x in oncokb.split(",")):
                    
                    if sift_score == "" or sift_pred in sift_tolerated:

                        if poly_pred == "" or poly_pred in polyphen_tolerated:
                            if not (poly_pred == "" and sift_pred == "" and clinvar == ""):
                                continue
            
            # add transcript to output
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

