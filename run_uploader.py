import os
import asyncio
from snip_warehouse import SnipUploader

DB_CONN_STRING = os.environ["SNIP_DB_CONN_STRING"]

async def main():
    snip_uploader = SnipUploader(DB_CONN_STRING)
    await snip_uploader.connect()
    f = open("genome_23andMe.txt", "r")
    await snip_uploader.upload(f, 1)


if __name__ == "__main__": 
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
    