import requests.adapters
from . import helpers
import requests
import json

def request_annotation(req_ls, proxies, verify, fields=None, token=None, json_out=None):
    """
    Request annotation for list of variants from genome nexus. For a list of possible fields to return see the genome nexus API.
    Some fields, like oncoKB, require a token for use and will not be returned if none is provided.
    """

    its = int(len(req_ls) / 4999) + 1
    print(f"{helpers.nice_time()} : Requesting annotation from Genome Nexus...")
    print(f"{helpers.nice_time()} : Making {its} requests...")

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
    
    resp_ls = []

    sesh = requests.Session()

    for i in range(its):
        last = min((i+1)*4999 - 1, len(req_ls) - 1)
        first = i*4999
        if first == last:
            ids = req_ls[first]
        else:
            ids = req_ls[first:(last+1)]
        
        retries = requests.adapters.Retry(total=5,
                                          backoff_factor=0.1)
        sesh.mount("https://", requests.adapters.HTTPAdapter(max_retries=retries))
        req_res = sesh.post(req_url, json=ids, headers=req_headers, params=payload, proxies=proxies, verify=verify)
        
        resp_status = req_res.status_code
        print(f"{helpers.nice_time()} : Request {i+1} returned exit code {resp_status}")
        resp = req_res.json()
        resp_ls.extend(resp)

    # old
    # req_res = requests.post(req_url, json=req_ls, headers=req_headers, params=payload, proxies=proxies, verify=verify)

    #resp_status = req_res.status_code
    #print(f"{helpers.nice_time()} : Request returned exit code {resp_status}")
    #resp_ls_dics = req_res.json()
    
    # for legacy reasons
    resp_ls_dics = resp_ls

    if json_out:
        with open(json_out, "w") as outfile:
            json.dump(resp_ls_dics, outfile, indent=4)
    
    return(resp_ls_dics)


def request_variant_info(req_ls, proxies, verify, max_gnomad):
    """
    For input variants, obtain frequency within the population from gnomAD, using the myvariant.info API.
    myvariant.info only accepts 1000 variants per call. A result is returned if the variant frequency is less than the specified maximum or there is no information.
    """

    its = int(len(req_ls) / 999) + 1
    
    print(f"{helpers.nice_time()} : Requesting gnomad info from myvariant.info...")

    print(f"{helpers.nice_time()} : Making {its} requests...")

    # myvariant.info expects chr<hgvsg> as variant format
    full_ids = [f"chr{x}" for x in req_ls]

    params = {"fields": "gnomad_genome"}

    req_url = "https://myvariant.info/v1/variant"

    req_headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
    
    resp_ls = []
    sesh = requests.Session()
    # request info in blocks of 999
    for i in range(its):
        last = min((i+1)*999 - 1, len(full_ids) - 1)
        first = i*999
        if first == last:
            payload = {"ids": full_ids[first]}
        else:
            payload = {"ids": full_ids[first:last]}
        success = False
        repeat = 0
        while not success and repeat < 5:
            retries = requests.adapters.Retry(total=5,
                                            backoff_factor=0.1)
            sesh.mount("https://", requests.adapters.HTTPAdapter(max_retries=retries))
            req_res = sesh.post(req_url, json=payload, headers=req_headers, params=params, proxies=proxies, verify=verify)
            
            resp_status = req_res.status_code
            print(f"{helpers.nice_time()} : Request {i+1} returned exit code {resp_status}")
            if resp_status == 200:
                success = True
            else:
                repeat += 1
        resp = req_res.json()
        resp_ls.extend(resp)


    # filter by population frequency
    filtered_resp_ls = []
    for item in resp_ls:
        try:
            if item["gnomad_genome"]["af"]["af"] <= max_gnomad:
                filtered_resp_ls.append(item)
        except KeyError:
            filtered_resp_ls.append(item)

    # create dictionary with hgvsg as keys
    indexed_resp = dict((item["query"][3:], item) for item in filtered_resp_ls)

    return(indexed_resp)
