import os
from snip_warehouse import SnipLoader
from snip_warehouse.schema import init_db

DB_CONN_STRING = os.environ["SNIP_DB_CONN_STRING"]
DB_NAME = 'snp-master'

init_db(DB_NAME)
snip_loader = SnipLoader(DB_CONN_STRING)
chr_suffixes = [i for i in range(3, 23)]
chr_suffixes += ["X", "Y", "MT"]
for chr_suffix in chr_suffixes:
    fname = f"refsnp-chr{chr_suffix}.json.bz2"
    print(f'Starting download for {fname}')

    if not os.path.isfile(fname):
        snip_loader.download_dbsnp_file(f"refsnp-chr{chr_suffix}.json.bz2",
                                        chr_suffix)
    else:
        print(f"refsnp-chr{chr_suffix}.json.bz2 previously downloaded - moving to load.")
    
    snip_loader.load_ref_snps(
        f"refsnp-chr{chr_suffix}.json.bz2", str(chr_suffix))
    
    os.system(f"rm refsnp-chr{chr_suffix}.json.bz2")

    # should_continue = input('Do you want to continue to the next chromosome? (Y/N)')
    # if should_continue.upper() == 'Y':
    #     break
