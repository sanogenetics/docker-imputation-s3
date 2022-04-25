from typing import Generator, Tuple
import re
from puretabix import VCFAccumulator, get_vcf_fsm
from puretabix.bgzip import BlockGZipReader
from puretabix.vcf import LINE_START, VCFLine


def generate_lines(filename):
    with open(filename, "rb") as infile:
        reader = BlockGZipReader(infile)
        infile.seek(0)

        accumulator = VCFAccumulator()
        vcf_fsm = get_vcf_fsm()
        for line in reader.generate_lines():

            try:
                vcf_fsm.run(line.decode() + "\n", LINE_START, accumulator)
                vcfline = accumulator.to_vcfline()
                accumulator.reset()
                if vcfline.comment_raw or vcfline.comment_key:
                    continue
                yield vcfline
            except Exception as e:
                print(filename)
                print(line)
                raise e


chr_order = {
    "chr1": 1,
    "chr2": 2,
    "chr3": 3,
    "chr4": 4,
    "chr5": 5,
    "chr6": 6,
    "chr7": 7,
    "chr8": 8,
    "chr9": 9,
    "chr10": 10,
    "chr11": 11,
    "chr12": 12,
    "chr13": 13,
    "chr14": 14,
    "chr15": 15,
    "chr16": 16,
    "chr17": 17,
    "chr18": 18,
    "chr19": 19,
    "chr20": 20,
    "chr21": 21,
    "chr22": 22,
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5": 5,
    "6": 6,
    "7": 7,
    "8": 8,
    "9": 9,
    "10": 10,
    "11": 11,
    "12": 12,
    "13": 13,
    "14": 14,
    "15": 15,
    "16": 16,
    "17": 17,
    "18": 18,
    "19": 19,
    "20": 20,
    "21": 21,
    "22": 22,
}


def get_matches(
    origs: Generator[VCFLine, None, None], subs: Generator[VCFLine, None, None]
) -> Generator[Tuple[VCFLine, VCFLine], None, None]:
    orig = next(origs, None)
    sub = next(subs, None)
    while orig and sub:
        # if they match by position
        orig_chrom_pos = (chr_order[orig.chrom], orig.pos)
        sub_chrom_pos = (chr_order[sub.chrom],sub.pos)

#        if orig.chrom != "chr1":
#            print(orig)
#            print(sub)
#            print(orig_chrom_pos)
#            print(sub_chrom_pos)
#            print(orig_chrom_pos == sub_chrom_pos)
#            print(orig_chrom_pos < sub_chrom_pos)
#            print(orig_chrom_pos > sub_chrom_pos)
#            1/0

        if orig_chrom_pos == sub_chrom_pos:
            yield orig, sub
            # move both along
            orig = next(origs, None)
            sub = next(subs, None)
        # if orig is before
        elif orig_chrom_pos < sub_chrom_pos:
            orig = next(origs, None)
        # if sub is first
        elif orig_chrom_pos > sub_chrom_pos:
            sub = next(subs, None)
        else:
            print(orig)
            print(sub)
            raise RuntimeError("Unable to move on")


def doublesplit(input:str, split1:str, split2:str):
    for i in input.split(split1):
        for j in i.split(split2):
            if j:
                yield j

def get_allele(line:VCFLine) -> str:
    options = [line.ref]+list(line.alt)

    gt =line.sample[0]["GT"]

    return "".join(sorted(options[int(i)] for i in doublesplit(gt, "/", "|")))

def check_if_score(orig:VCFLine, sub:VCFLine):
    # orig and sub match to same chrom ans possition
    assert chr_order[orig.chrom] == chr_order[sub.chrom]
    assert orig.pos == sub.pos

    # check ref matches
    assert orig.ref == sub.ref
    # check alt matches
    #assert orig.alt == sub.alt

    # check if the sub has an imputed flag
    if "IMP" in sub.info:
        # calculate what alleles are spcified
        origallele = get_allele(orig)
        suballele = get_allele(sub)
        return origallele == suballele
    else:
        return None


if __name__ == "__main__":

    # calculate a cross validation score
    hit = 0
    miss = 0

    # for each imputed subset
    printed = set()
    for subfilename in (
        "/home/adam/bioinf/beagle/out/out.p01.vcf.gz",
        #"/home/adam/bioinf/beagle/out/out.p02.vcf.gz",
        #"/home/adam/bioinf/beagle/out/out.p03.vcf.gz",
        #"/home/adam/bioinf/beagle/out/out.p04.vcf.gz",
        #"/home/adam/bioinf/beagle/out/out.p05.vcf.gz",
    ):
        # load the original
        origs = generate_lines("/home/adam/bioinf/beagle/in/sano-b38.vcf.gz")
        subs = generate_lines(subfilename)
        # look for variants that are in original and imputed
        for orig, sub in get_matches(origs, subs):
            if orig.chrom not in printed:
                print(orig)
                print(sub)
                printed.add(orig.chrom)

            try:
                match = check_if_score(orig, sub)
            except Exception as e:
                print(orig)
                print(sub)
                raise e

            if match == None:
                pass
            elif match:
                hit += 1
            else:
                miss += 1

    # count how many are correct / incorrect
    print(hit, miss, hit/(hit+miss))
