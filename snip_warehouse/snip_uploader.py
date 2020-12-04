import asyncpg
import csv


class SnipUploader:
    def __init__(self, db_conn_string):
        self.db_conn_string = db_conn_string

    async def connect(self):
        self.pool = await asyncpg.create_pool(self.db_conn_string)
        self.conn = await self.pool.acquire()

    async def upload(self, iterable_buff, user_id):
        """ `iterable_buff`: Iterable
             -  will never be loaded entirely into memory!
             `user_id`: int
             -  the id of the user uploading SNP data
        """
        await self._upload(iterable_buff, user_id)

    async def _upload(self, buff, user_id=1):
        cr = csv.reader(filter(lambda ln: not ln.startswith("#"), buff),
                        delimiter="\t")
        next(cr)
        record_gen = ((int(r[0].replace("rs", "").replace("i", "")),
                       # int(r[1]),
                       r[3],
                       user_id,) for r in cr)
        await self.conn.copy_records_to_table(
            "user_ref_snps",
            columns=("ref_snp_id", "genotype", "user_id",),
            records=record_gen)
