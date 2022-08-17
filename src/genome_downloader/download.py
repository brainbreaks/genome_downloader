import tarfile
import requests
import httplib2
import argparse
import zipfile
import tqdm
import tempfile
import re
import shutil
import os
import glob
import gzip
from pathlib import Path
import subprocess

try:
    from genome_downloader.gff_longest_transcript import find_longest_transcript
    is_longest_trancript_generated = True
except ImportError:
    is_longest_trancript_generated = False
    print("gff_longest_transcript.py not found. Bed file with longest transcript per each gene will not be created...")


def url_exist(url):
    h = httplib2.Http()
    resp = h.request(url, 'HEAD')
    return int(resp[0]['status']) < 400

def download_ucsc_table(genome, track, table, dest=None, description=None, overwrite=False):
    data = {"jsh_pageVertPos": 0, "clade": "mammal", "org": "Mouse", "db": genome, "hgta_group": "varRep",
        "hgta_track": track, "hgta_table": table, "hgta_regionType": "genome", "position": "chr12:56694976-56714605",
        "hgta_outputType": "primaryTable", "boolshad.sendToGalaxy": 0, "boolshad.sendToGreat": 0,
        "hgta_outFileName": "output.tsv", "hgta_compressType": "none", "hgta_doTopSubmit": "get output"
    }

    path = download_file("https://genome.ucsc.edu/cgi-bin/hgTables", data=data, dest=dest, description=description, overwrite=overwrite)

    return path

def download_file(url, headers=None, data=None, dest=None, description=None, overwrite=False, compressed=False):
    if description is None:
        description = 'Downloading file "{}" ==> "{}"...'.format(url, dest)

    if not overwrite and dest and os.path.isfile(dest):
        print(description)
        return dest

    if re.match("^(http|ftp)", url):
        response = None
        if data is None:
            response = requests.get(url, stream=True)
        else:
            response = requests.post(url, data=data, stream=True)
        total_size_in_bytes = int(response.headers.get('content-length', 0))
        # size = response.headers.get('content-length', 0)
        # MBFACTOR = float(1 << 20)
        # print('{:<40}: {:.2f} MB'.format('FILE SIZE', int(size) / MBFACTOR))

        progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='B', unit_scale=True, unit_divisor=1024)
        with tempfile.NamedTemporaryFile(mode="wb", delete=False) as file:
            path = file.name
            progress_bar.set_description(description)
            for data in response.iter_content(1024):
                progress_bar.update(len(data))
                file.write(data)
            progress_bar.close()
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            print("ERROR, something went wrong")
    else: # local file
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        shutil.copy2(url, tmp.name)
        path = tmp.name

    if compressed:
        with gzip.open(path, 'rb') as f_in:
            with tempfile.NamedTemporaryFile(mode="wb", delete=False) as f_out:
                shutil.copyfileobj(f_in, f_out)
                path = f_out.name

    if dest:
        dest = os.path.abspath(dest)
        if not os.path.exists(os.path.dirname(dest)):
            os.makedirs(os.path.dirname(dest))

        shutil.copy2(path, dest)
        path = dest

    return path



def download_raw_genome(url, dest, overwrite=False):
    if not overwrite and dest and os.path.isfile(dest):
        print('Downloading raw genome "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
        return

    if not os.path.exists(os.path.dirname(dest)):
        os.makedirs(os.path.dirname(dest))

    path = download_file(url, description='Downloading raw genome "{}"\n'.format(url))

    dest_chromFa = "{}_{}".format(dest, "chromFa")
    with tarfile.open(path, mode="r:gz") as tf:
        progress_bar = tqdm.tqdm(tf.getmembers(), desc="Extracting raw genome archive contents into '{}'\n".format(dest_chromFa))
        for member in progress_bar:
            try:
                tf.extract(member, dest_chromFa)
            except tarfile.TarError as e:
                pass

    print("Joining chromosomes into single fasta file '{}'\n".format(dest))
    with open(dest, 'w') as outfile:
        # Iterate through list
        for path in Path(dest_chromFa).rglob('*.fa'):
            with open(str(path)) as infile:
                outfile.write(infile.read())
            outfile.write("\n")

    shutil.rmtree(dest_chromFa)



def download_bowtie2_index(url, dest, genome, overwrite=False):
    if not overwrite and dest and len(glob.glob(os.path.join(dest, "{}.*.bt2".format(genome)))) > 0:
        print('Downloading bowtie2 index "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
        return

    path = download_file(url, description="Downloading bowtie index '{}'\n".format(url))

    with zipfile.ZipFile(path, 'r') as zf:
        for member in tqdm.tqdm(zf.infolist(), desc="Extracting bowtie2 index into '{}'\n".format(dest)):
            try:
                zf.extract(member, dest)
            except zipfile.error as e:
                pass

def build_fasta_index(path):
    print('Building fasta index for "{}" ==> "{}i"...'.format(path, path))
    "samtools faidx hg38.fa"
    subprocess.run(["samtools", "faidx", path])

def download_genome(genome, path):
    download_ucsc_table(genome=genome, table="rmsk", track="rmsk", dest=os.path.join(path, "{genome}/annotation/ucsc_repeatmasker.tsv".format(genome=genome)))

    download_file("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/chromInfo.txt.gz".format(genome=genome), dest=os.path.join(path, "{genome}/annotation/ChromInfo.txt".format(genome=genome)), compressed=True)
    download_file("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/cytoBand.txt.gz".format(genome=genome), dest=os.path.join(path, "{genome}/annotation/cytoBand.txt".format(genome=genome)), compressed=True)
    if url_exist("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/chromFa.tar.gz".format(genome=genome)):
        download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/chromFa.tar.gz".format(genome=genome), dest=os.path.join(path, "{genome}/{genome}.fa".format(genome=genome)))
    elif url_exist("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chromFa.tar.gz".format(genome=genome)):
        download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chromFa.tar.gz".format(genome=genome), dest=os.path.join(path, "{genome}/{genome}.fa".format(genome=genome)))
    else:
        raise FileExistsError("chromFa.tar.gz url couldn't be accessed")

    build_fasta_index(os.path.join(path, "{genome}/{genome}.fa").format(genome=genome))

    download_file("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes".format(genome=genome), dest=os.path.join(path, "{genome}/annotation/{genome}.chrom.sizes".format(genome=genome)))
    download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/{genome}.zip".format(genome=genome), genome=genome, dest=os.path.join(path, genome))

    download_file("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.ncbiRefSeq.gtf.gz".format(genome=genome), dest=os.path.join(path, "{genome}/annotation/{genome}.ncbiRefSeq.gtf.gz".format(genome=genome)))
    download_file("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.refGene.gtf.gz".format(genome=genome), dest=os.path.join(path, "{genome}/annotation/{genome}.refGene.gtf.gz".format(genome=genome)))
    print("Creating annotation file from {gtf}...".format(gtf=os.path.join(path, "{genome}/annotation/{genome}.refGene.gtf.gz".format(genome=genome))))

    if is_longest_trancript_generated:
        find_longest_transcript(os.path.join(path, "{genome}/annotation/{genome}.refGene.gtf.gz".format(genome=genome)), os.path.join(path, "{genome}/annotation/refGene.bed".format(genome=genome)), clip_start=50, clip_strand_specific=True)

def cli():
    parser = argparse.ArgumentParser(description='Download groseq dependencies')
    parser.add_argument('data', choices=['mm9', 'mm10', 'hg19'], help="""Should be one of the following: \n
    mm9 - mm9 model files\n  
    mm10 - mm10 model files\n
    hg19 - hg19 model files
    """)
    parser.add_argument('path', nargs="?", default=".", help="Path where the files will be downloaded (default: current directory)")

    args = parser.parse_args()

    if not os.path.exists(args.path):
        os.makedirs(args.path)

    download_genome(args.data, args.path)

if __name__ == '__main__':
    cli()