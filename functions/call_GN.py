from . import helpers
import requests
import json

def request_annotation(req_ls, proxies, verify, fields=None, token=None, json_out=None):



    print(f"{helpers.nice_time()} : Requesting annotation from Genome Nexus...")

    payload = {}
    if token:
        payload["token"] = json.dumps(token)
        #payload = {"token": json.dumps(token),
         #          "fields": fields}
    if fields:
        payload["fields"] = fields

    req_url = f"https://www.genomenexus.org/annotation"

    req_headers = {"Content-Type": "application/json",
                "Accept": "application/json"}
    
    req_res = requests.post(req_url, json=req_ls, headers=req_headers, params=payload, proxies=proxies, verify=verify)

    resp_status = req_res.status_code
    print(f"{helpers.nice_time()} : Request returned exit code {resp_status}")
    resp_ls_dics = req_res.json()

    if json_out:
        with open(json_out, "w") as outfile:
            json.dump(resp_ls_dics, outfile, indent=4)
    
    return(resp_ls_dics)


def request_variant_info(req_ls, proxies, verify, max_gnomad):

    its = int(len(req_ls) / 999) + 1
    
    print(f"{helpers.nice_time()} : Requesting gnomad info from myvariant.info...")

    print(f"{helpers.nice_time()} : Making {its} requests...")

    full_ids = [f"chr{x}" for x in req_ls]

    params = {"fields": "gnomad_genome"}

    req_url = "https://myvariant.info/v1/variant"

    req_headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
    
    resp_ls = []

    for i in range(its):
        last = min((i+1)*999 - 1, len(full_ids) - 1)
        first = i*999
        if first == last:
            payload = {"ids": full_ids[first]}
        else:
            payload = {"ids": full_ids[first:last]}
        req_res = requests.post(req_url, json=payload, headers=req_headers, params=params, proxies=proxies, verify=verify)
        
        resp_status = req_res.status_code
        print(f"{helpers.nice_time()} : Request {i+1} returned exit code {resp_status}")
        resp = req_res.json()
        resp_ls.extend(resp)


    filtered_resp_ls = []
    for item in resp_ls:
        try:
            if item["gnomad_genome"]["af"]["af"] <= max_gnomad:
                filtered_resp_ls.append(item)
        except KeyError:
            continue
            
    indexed_resp = dict((item["query"][3:], item) for item in filtered_resp_ls)



    return(indexed_resp)
