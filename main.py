from aiohttp import web
import asyncpg
import ujson as json


async def get_user_diseases(request):
    # TODO: Pull this off of Auth context variable
    user_id = int(request.match_info["user_id"])
    # TODO: Store citation_csv as an Array
    async with app['conn_pool'].acquire() as conn:
        significance_group_future = conn.fetch("""
            SELECT array_agg(disease_name_csv) disease_name_ls,
                   array_agg(array_to_string(citation_list, ',')) citation_ls,
                   array_agg(s.genotype) genotype_ls,
                   array_agg(a.ins_seq) ins_seq_ls,
                   array_agg(s.ref_snp_id) ref_snp_id_ls,
                   Count(*) count,
                   clinical_significance_csv
            FROM ref_snp_allele_clin_diseases d
            INNER JOIN ref_snp_alleles a
            ON a.id = d.ref_snp_allele_id
            INNER JOIN user_ref_snps s ON s.ref_snp_id = a.ref_snp_id
            WHERE s.user_id = $1
             AND ((s.genotype iLIKE '%' || a.ins_seq  || '%'
                  AND a.ins_seq != ''
                  AND a.del_seq != '')
                  OR (s.genotype iLIKE '%I%'
                      AND a.ins_seq != ''
                      AND a.del_seq = '')
                  OR (s.genotype iLIKE '%D%'
                      AND a.ins_seq = ''
                      AND a.del_seq != ''))
            GROUP BY clinical_significance_csv""", user_id)
        clin_sig_groups = {
            row['clinical_significance_csv']: [
                {
                   "disease_name": d,
                   "ins_seq": alt,
                   "genotype": geno,
                   "citation_list": cit.split(","),
                   "ref_snp_id": rsnp_id
                } for d, cit, alt, geno, rsnp_id in zip(
                    row['disease_name_ls'],
                    row['citation_ls'],
                    row['ins_seq_ls'],
                    row['genotype_ls'],
                    row['ref_snp_id_ls']
                )]
            for row in await significance_group_future}
        res = web.Response(body=json.dumps(clin_sig_groups))
        res.headers['Content-Type'] = "application/json"
        return res


def init_routes(app):
    resource = app.router.add_resource("/users/{user_id}/diseases")
    resource.add_route("GET", get_user_diseases)


async def init_conn_pools(app):
    num_bulk_upload_connections = app['config']['num_bulk_upload_connections']
    num_fetch_connections = app['config']['num_connections']
    bulk_upload_conn_pool = await asyncpg.create_pool(
        user="SeanH",
        database="snvs_chr_1",
        max_size=num_bulk_upload_connections,
        loop=app.loop)
    # TODO: Set connection kwargs to ignore integrity checks?
    conn_pool = await asyncpg.create_pool(
        user="SeanH",
        database="snvs_chr_1",
        max_size=num_fetch_connections,
        loop=app.loop)
    app['bulk_upload_conn_pool'] = bulk_upload_conn_pool
    app['conn_pool'] = conn_pool


async def close_conn_pools(app):
    await app['bulk_upload_conn_pool'].close()
    await app['conn_pool'].close()


app = web.Application()
app['config'] = {
    "num_connections": 15,
    "num_bulk_upload_connections": 10
}
init_routes(app)
app.on_startup.append(init_conn_pools)
app.on_cleanup.append(close_conn_pools)
web.run_app(app, host="127.0.0.1", port=8080)
