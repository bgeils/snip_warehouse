import asyncpg
import asyncio
import ftplib
import bz2       
from multiprocessing import cpu_count, Pool
import time
import threading
from typing import List
import ujson as json

from .types import (
    RefSnpCopyFromData,
    RefSnpAllele,
    RefSnpAlleleFreqStudy,
    RefSnpAlleleClinDisease)

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
'''
This class includes functions to download, parse and load NIH SNP files.

Each file downloaded corresponds to a chromosome ()
'''
class SnipLoader:
    def __init__(self, db_conn_string):
        self.db_conn_string = db_conn_string
        self.file_blocksize = 1024 * 1024

    def download_dbsnp_file(self, dbsnp_filename, chromosome):
        self.chromosome = str(chromosome)
        self.dbsnp_filename = dbsnp_filename
        with open(dbsnp_filename, "wb") as fp:
            ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
            ftp.login()
            ftp.cwd("snp/.redesign/latest_release/JSON")
            file_size = ftp.size(dbsnp_filename)
            size_gb = round(file_size / (1024**3), 2)
            size_mb = int(file_size / (1024**2))
            print(f"Filesize: {size_gb} GB")
            sock = None
            while not sock:  # Try to open the socket conn w/ server
                print("Trying to establish FTP conn")
                sock = ftp.transfercmd(f"RETR {dbsnp_filename}")
                time.sleep(5)

            def download():
                transferred, blocks = 0, 0
                while True:
                    byte_chunk = sock.recv(1024*1024*8)
                    if not byte_chunk:
                        break
                    blocks += 1
                    transferred += len(byte_chunk)
                    transferred_mb = round(transferred / 1024 / 1024, 2)
                    if blocks % 1000 == 0:
                        printProgressBar(int(transferred_mb), size_mb, prefix = 'Progress:', suffix = 'Complete', length = 50, decimals = 3)
                    fp.write(byte_chunk)
                sock.close()
                fp.close()
                print('\n')
            t = threading.Thread(target=download)
            t.start()
            while t.is_alive():
                t.join(60)
                ftp.voidcmd("NOOP")

    def load_ref_snps(self, dbsnp_filename, chromosome):
        self.chromosome = chromosome
        num_processes = cpu_count()
        print(f"Found '{num_processes}' CPUs")
        loop = asyncio.get_event_loop()
        with Pool(num_processes) as pool:
            with bz2.open(dbsnp_filename) as bz2_fp:
                print("Mapping...")
                copy_from_data_iter = pool.imap_unordered(
                    self._generate_parsed_data,
                    bz2_fp,
                    1024)
                loop.run_until_complete(self._load(copy_from_data_iter))

    async def _load(self, parsed_data_iter):
        conn = await asyncpg.connect(self.db_conn_string)
        table_names = [table_name for table_name in RefSnpCopyFromData._fields]
        row_buff_dict = {table_name: [] for table_name in table_names}
        buff_size = 0
        for copy_from_data in parsed_data_iter:
            for table_name in table_names:
                for copy_from_row in getattr(copy_from_data, table_name):
                    row_buff_dict[table_name].append(copy_from_row)
            buff_size += 1
            if buff_size % 5000 == 0:  # Dump
                print(f"Dumping (SNPs Processed: {buff_size})")
                await self._dump_buffer(row_buff_dict, conn)
                row_buff_dict = {table_name: [] for table_name in table_names}
                print("Done.")
        await conn.close()

    async def _dump_buffer(self, row_buff_dict, conn):
        for table_name in row_buff_dict.keys():
            records = row_buff_dict[table_name]
            if not records:
                continue
            try: 
                await conn.copy_records_to_table(
                    table_name,
                    records=records,
                    columns=records[0]._fields)
            except asyncpg.exceptions.UniqueViolationError as e:
                continue

    def _generate_parsed_data(self, raw_line) -> RefSnpCopyFromData:
        rsnp_json = json.loads(raw_line)
        ref_snp_id = int(rsnp_json['refsnp_id'])
        rsnp_placements = rsnp_json['primary_snapshot_data'][
                                'placements_with_allele']
        copy_from_data = RefSnpCopyFromData(
            ref_snp_alleles=[],
            ref_snp_allele_freq_studies=[],
            ref_snp_allele_clin_diseases=[])
        if not rsnp_placements:
            return copy_from_data
        allele_data = self._find_alleles_from_assembly(rsnp_placements)
        if not allele_data:
            return copy_from_data
        variant_ref_snp_alleles = self._get_variant_alleles(self.chromosome,
                                                            allele_data,
                                                            ref_snp_id)
        if not variant_ref_snp_alleles:
            return copy_from_data
        for allele in variant_ref_snp_alleles:
            allele_idx = allele.ref_snp_allele_idx
            allele_annotation = rsnp_json['primary_snapshot_data'][
                'allele_annotations'][allele_idx]
            gene_locii = self._parse_gene_locii(allele_annotation)
            freq_studies = self._parse_freq_studies(allele_annotation,
                                                    ref_snp_id, allele_idx)
            clin_diseases = self._parse_clin_diseases(allele_annotation,
                                                      ref_snp_id, allele_idx)
            copy_from_data = self._update_copy_from_data(
               copy_from_data, allele, freq_studies, clin_diseases, gene_locii)
        return copy_from_data

    @staticmethod
    def _find_alleles_from_assembly(rsnp_placements,
                                    assembly_name="GRCh38"):
        for rsnp_placement in rsnp_placements:
            annot = rsnp_placement.get('placement_annot')
            if not annot or not annot.get('seq_id_traits_by_assembly'):
                return
            assembly_info_ls = annot['seq_id_traits_by_assembly']
            assembly_info = assembly_info_ls[0]
            this_assembly_name = assembly_info.get("assembly_name")
            if assembly_name in this_assembly_name:
                return rsnp_placement['alleles']

    @staticmethod
    def _get_variant_alleles(chromosome, alleles,
                             ref_snp_id) -> List[RefSnpAllele]:
        variant_alleles = []
        for idx, allele in enumerate(alleles):
            spdi = allele['allele']['spdi']
            ins, delete = spdi['inserted_sequence'], spdi['deleted_sequence']
            if ins != delete:
                variant_alleles.append(RefSnpAllele(
                    ins_seq=ins,
                    del_seq=delete,
                    position=spdi['position'],
                    ref_snp_allele_idx=idx,
                    chromosome=chromosome,
                    ref_snp_id=ref_snp_id,
                    gene_locii=None))
        return variant_alleles

    @staticmethod
    def _update_copy_from_data(copy_from_data: RefSnpCopyFromData,
                               allele: RefSnpAllele,
                               freq_studies: List[RefSnpAlleleFreqStudy],
                               clin_diseases: List[RefSnpAlleleClinDisease],
                               gene_locii: List[str]):
        copy_from_data.ref_snp_alleles.append(
            RefSnpAllele(
                del_seq=allele.del_seq,
                ins_seq=allele.ins_seq,
                position=allele.position,
                ref_snp_allele_idx=allele.ref_snp_allele_idx,
                chromosome=allele.chromosome,
                ref_snp_id=allele.ref_snp_id,
                gene_locii=gene_locii))
        copy_from_data.ref_snp_allele_freq_studies.extend(freq_studies)
        copy_from_data.ref_snp_allele_clin_diseases.extend(clin_diseases)
        return copy_from_data

    @staticmethod
    def _parse_freq_studies(allele_annotation, ref_snp_id, allele_idx):
        return [RefSnpAlleleFreqStudy(
            name=freq['study_name'],
            allele_count=freq['allele_count'],
            total_count=freq['total_count'],
            ref_snp_allele_idx=allele_idx,
            ref_snp_id=ref_snp_id)
            for freq in allele_annotation['frequency'] or []]

    @staticmethod
    def _parse_clin_diseases(allele_annotation, ref_snp_id, allele_idx):
        return [RefSnpAlleleClinDisease(
            citation_list=clin['citations'],
            disease_name_csv=",".join(clin['disease_names']),
            clinical_significance_csv=",".join(clin['clinical_significances']),
            ref_snp_allele_idx=allele_idx,
            ref_snp_id=ref_snp_id)
                for clin in allele_annotation['clinical']]

    @staticmethod
    def _parse_gene_locii(allele_annotation):
        assembly_annotation = allele_annotation['assembly_annotation']
        return set(
            [gene['locus'] for gene in
                (assembly_annotation[0]['genes']
                 if assembly_annotation else [])])
